"""
organisms.py
Species / organism resolution.

THE HEADLINE BUGFIX
-------------------
v2.x resolved a species by fetching ``https://rest.kegg.jp/list/organism`` at
call time.  That had three independent failure modes:

  1. No offline data.  Behind a firewall (or when KEGG throttles, which it does
     aggressively) the fetch raised and *every* pathview() call died before it
     reached any pathway logic.  This is the "can't find the organism list link"
     symptom.
  2. Whole-column equality matching.  The organism list stores the name as
     ``Homo sapiens (human)``.  v2.x tested ``query in (p.lower() for p in
     parts)``, an exact match against a whole column, so "human", "Homo
     sapiens", "homo sapiens" and every other natural spelling failed even when
     the fetch succeeded.  Only the code itself ("hsa") ever matched.
  3. Wrong column index.  KEGG's organism list is
     ``ktax_id, code, name, lineage``; v2.x read ``parts[1]`` as the code, which
     is right for that endpoint but silently wrong for the ``/list/organism``
     variants that omit the T-number, returning a taxonomy id as a species code.

The fix: ship the organism table (all 10,718 KEGG organisms, the same table R
pathview distributes as ``korg``) inside the package and resolve against it in
memory.  Zero network calls, works behind a firewall, and matches on code,
T-number, taxonomy id, scientific name, common name, parenthesised alias,
genus abbreviation and case-insensitive substrings.  KEGG is consulted only to
*refresh* the bundled table, never to satisfy a lookup.

Public API
----------
  get_species_code / kegg_species_code : resolve a species -> SpeciesInfo
  list_organisms                       : the full table as a DataFrame
  search_organisms                     : fuzzy search
  refresh_organism_table               : re-download from KEGG (optional)
"""

from __future__ import annotations

import difflib
import gzip
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import polars as pl

from .constants import KEGG_BASE
from .errors import NetworkError, SpeciesNotFoundError

_DATA = Path(__file__).parent / "data"
_KORG = _DATA / "korg.tsv.gz"
_BODS = _DATA / "bods.tsv.gz"


# ---------------------------------------------------------------------------
# SpeciesInfo
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SpeciesInfo:
    """Immutable KEGG organism record."""

    kegg_code: str
    scientific_name: str = ""
    common_name: str = ""
    ktax_id: str = ""
    tax_id: str = ""
    entrez_gnodes: bool = True
    kegg_geneid: str | None = None
    ncbi_geneid: str | None = None
    ncbi_proteinid: str | None = None
    uniprot: str | None = None

    @property
    def display_name(self) -> str:
        if self.common_name and self.scientific_name:
            return f"{self.scientific_name} ({self.common_name})"
        return self.scientific_name or self.kegg_code

    def __str__(self) -> str:                              # pragma: no cover
        return f"{self.kegg_code}: {self.display_name}"


# Orthology pseudo-organism: KEGG "ko" maps are species-agnostic.
KO_SPECIES = SpeciesInfo(
    kegg_code="ko",
    scientific_name="KEGG Orthology",
    common_name="ortholog reference",
    entrez_gnodes=False,
    kegg_geneid="K01488",
)

# Curated vernacular aliases.
#
# Substring search alone is ambiguous for common names: "mouse" matches both
# *Mus musculus* ("house mouse") and *Microcebus murinus* ("gray mouse lemur"),
# and a naive first-hit wins the wrong one.  These are the organisms a user
# typing a bare vernacular name unambiguously means.
_ALIASES: dict[str, str] = {
    "mouse": "mmu", "mice": "mmu", "murine": "mmu",
    "rat": "rno", "human": "hsa", "man": "hsa",
    "fly": "dme", "fruit fly": "dme", "drosophila": "dme",
    "worm": "cel", "nematode": "cel", "c elegans": "cel", "c. elegans": "cel",
    "zebrafish": "dre", "zebra fish": "dre",
    "yeast": "sce", "budding yeast": "sce", "baker's yeast": "sce",
    "fission yeast": "spo",
    "arabidopsis": "ath", "thale cress": "ath",
    "e coli": "eco", "e. coli": "eco", "ecoli": "eco",
    "b subtilis": "bsu", "bacillus subtilis": "bsu",
    "chicken": "gga", "cow": "bta", "cattle": "bta", "bovine": "bta",
    "pig": "ssc", "swine": "ssc", "porcine": "ssc",
    "dog": "cfa", "canine": "cfa", "cat": "fca", "feline": "fca",
    "sheep": "oas", "horse": "eca", "equine": "eca",
    "rabbit": "ocu", "macaque": "mcc", "rhesus": "mcc",
    "chimp": "ptr", "chimpanzee": "ptr",
    "maize": "zma", "corn": "zma", "rice": "osa",
    "wheat": "taes", "soybean": "gmx", "tomato": "sly",
    "frog": "xtr", "xenopus": "xtr",
    "tb": "mtu", "m tuberculosis": "mtu", "mycobacterium tuberculosis": "mtu",
    "p aeruginosa": "pae", "s aureus": "sau", "staph": "sau",
    "salmonella": "sey", "malaria": "pfa", "plasmodium": "pfa",
}

# Enzyme / reaction / compound reference maps
_PSEUDO = {
    "ko": KO_SPECIES,
    "ec": SpeciesInfo(kegg_code="ec", scientific_name="Enzyme Commission",
                      common_name="enzyme reference", entrez_gnodes=False),
    "rn": SpeciesInfo(kegg_code="rn", scientific_name="KEGG Reaction",
                      common_name="reaction reference", entrez_gnodes=False),
    "map": SpeciesInfo(kegg_code="map", scientific_name="KEGG reference pathway",
                       common_name="reference", entrez_gnodes=False),
    "cpd": SpeciesInfo(kegg_code="cpd", scientific_name="KEGG Compound",
                       common_name="compound reference", entrez_gnodes=False),
}


# ---------------------------------------------------------------------------
# Table loading
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def _load_table() -> pl.DataFrame:
    """Load the bundled organism table (cached for the process lifetime)."""
    override = _user_table_path()
    src = override if override.exists() else _KORG
    with gzip.open(src, "rt", encoding="utf-8") as fh:
        df = pl.read_csv(
            fh.read().encode(),
            separator="\t",
            infer_schema_length=0,      # everything as str: ids are not numbers
        )
    return df.with_columns(
        pl.col("kegg_code").str.to_lowercase().alias("_code_lc"),
        pl.col("scientific_name").str.to_lowercase().alias("_sci_lc"),
        pl.col("common_name").str.to_lowercase().alias("_com_lc"),
    )


def _user_table_path() -> Path:
    from .cache import cache_dir
    return cache_dir() / "korg.refreshed.tsv.gz"


def list_organisms(
    pattern: str | None = None,
    with_genome_only: bool = False,
) -> pl.DataFrame:
    """
    Return the bundled KEGG organism table.

    Parameters
    ----------
    pattern:          Optional case-insensitive regex filtered against the KEGG
                      code, scientific name and common name.
    with_genome_only: Keep only organisms with Entrez gene nodes.

    Examples
    --------
    >>> list_organisms("mouse").select(["kegg_code", "scientific_name"])
    """
    df = _load_table().drop(["_code_lc", "_sci_lc", "_com_lc"])
    if with_genome_only:
        df = df.filter(pl.col("entrez_gnodes") == "1")
    if pattern:
        rx = f"(?i){pattern}"
        df = df.filter(
            pl.col("kegg_code").str.contains(rx)
            | pl.col("scientific_name").str.contains(rx)
            | pl.col("common_name").str.contains(rx)
        )
    return df


#: KEGG codes for the reference organisms named in :data:`_ALIASES`.
#: Used to rank search results towards the organism most people mean.
_REFERENCE_CODES: frozenset[str] = frozenset(_ALIASES.values())


def search_organisms(query: str, n: int = 10,
                     limit: int | None = None) -> pl.DataFrame:
    """
    Fuzzy-search the organism table; useful after a SpeciesNotFoundError.

    Results are ranked, not merely filtered.  A plain substring filter puts
    *Yersinia enterocolitica* above *Escherichia coli* for the query "coli",
    because "enterocolitica" contains "coli" too.  Scoring therefore rewards
    exact fields, then whole-word matches, then prefixes, and prefers short
    names so reference strains outrank long strain designations.

    Returns a DataFrame with a ``score`` column, highest first.
    """
    if limit is not None:
        n = int(limit)
    n = max(1, int(n))
    df = _load_table()
    q = str(query or "").lower().strip()
    if not q:
        return df.head(0).drop(["_code_lc", "_sci_lc", "_com_lc"])

    word = rf"(?:^|[^a-z0-9]){re.escape(q)}(?:[^a-z0-9]|$)"

    # Match score first: only rows that actually match the query survive.
    # The reference-organism bonus must not be part of this, or every
    # reference organism would "match" any query at all.
    match = (
        pl.when(pl.col("_code_lc") == q).then(1000).otherwise(0)
        + pl.when(pl.col("_sci_lc") == q).then(900).otherwise(0)
        + pl.when(pl.col("_com_lc") == q).then(850).otherwise(0)
        + pl.when(pl.col("_sci_lc").str.contains(word)).then(400).otherwise(0)
        + pl.when(pl.col("_com_lc").str.contains(word)).then(380).otherwise(0)
        + pl.when(pl.col("_sci_lc").str.starts_with(q)).then(300).otherwise(0)
        + pl.when(pl.col("_code_lc").str.starts_with(q)).then(200).otherwise(0)
        + pl.when(pl.col("_sci_lc").str.contains(re.escape(q))).then(100).otherwise(0)
        + pl.when(pl.col("_com_lc").str.contains(re.escape(q))).then(90).otherwise(0)
    ).alias("_match")

    scored = df.with_columns(match).filter(pl.col("_match") > 0).with_columns(
        (
            pl.col("_match")
            # Reference organisms (those the curated alias table points at)
            # outrank strain-specific entries: "coli" wants E. coli K-12, not
            # an arbitrary sequenced isolate.
            + pl.when(pl.col("kegg_code").is_in(_REFERENCE_CODES)).then(150).otherwise(0)
            # Shorter names are usually the reference entry; break ties there.
            - (pl.col("_sci_lc").str.len_chars() / 10).cast(pl.Int32)
        ).alias("score")
    ).drop("_match")

    if scored.is_empty():                      # last resort: edit distance
        close = difflib.get_close_matches(q, df["_sci_lc"].to_list(), n=n, cutoff=0.6)
        if not close:
            # Try the first token: "Homo sapein" -> "homo" still finds hsa.
            head = q.split()[0] if q.split() else q
            close = difflib.get_close_matches(head, df["_sci_lc"].to_list(),
                                              n=n, cutoff=0.5)
        if not close:
            return (df.head(0).with_columns(pl.lit(0).alias("score"))
                      .drop(["_code_lc", "_sci_lc", "_com_lc"]))
        scored = (df.filter(pl.col("_sci_lc").is_in(close))
                    .with_columns(
                        (pl.lit(50)
                         + pl.when(pl.col("kegg_code").is_in(_REFERENCE_CODES))
                           .then(20).otherwise(0)).alias("score")))

    return (scored.sort(["score", "kegg_code"], descending=[True, False])
                  .head(n)
                  .drop(["_code_lc", "_sci_lc", "_com_lc"]))


# ---------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------

def _row_to_info(row: dict) -> SpeciesInfo:
    def clean(v: object) -> str | None:
        s = str(v).strip() if v is not None else ""
        return s or None

    return SpeciesInfo(
        kegg_code=str(row["kegg_code"]).strip(),
        scientific_name=str(row.get("scientific_name") or "").strip(),
        common_name=str(row.get("common_name") or "").strip(),
        ktax_id=str(row.get("ktax_id") or "").strip(),
        tax_id=str(row.get("tax_id") or "").strip(),
        entrez_gnodes=str(row.get("entrez_gnodes") or "1").strip() in ("1", "True", "true"),
        kegg_geneid=clean(row.get("kegg_geneid")),
        ncbi_geneid=clean(row.get("ncbi_geneid")),
        ncbi_proteinid=clean(row.get("ncbi_proteinid")),
        uniprot=clean(row.get("uniprot")),
    )


@lru_cache(maxsize=512)
def get_species_code(species: str = "hsa") -> SpeciesInfo:
    """
    Resolve *species* to a :class:`SpeciesInfo`.

    Accepts a KEGG code (``hsa``), T-number (``T01001``), NCBI taxonomy id
    (``9606``), scientific name (``Homo sapiens``), common name (``human``),
    or a name with a parenthesised alias (``Homo sapiens (human)``).  Matching
    is case-insensitive and whitespace-tolerant.

    Resolution is entirely offline: it reads the organism table bundled with
    the package.  No network access is attempted at any point.

    Raises
    ------
    SpeciesNotFoundError
        With near-miss suggestions drawn from the table.

    Examples
    --------
    >>> get_species_code("human").kegg_code
    'hsa'
    >>> get_species_code("Mus musculus").kegg_code
    'mmu'
    >>> get_species_code("9606").kegg_code
    'hsa'
    """
    if species is None:
        raise SpeciesNotFoundError("None")

    raw = str(species).strip()
    if not raw:
        raise SpeciesNotFoundError("")

    key = re.sub(r"\s+", " ", raw.lower()).strip()
    if key in _PSEUDO:
        return _PSEUDO[key]

    df = _load_table()

    # 0. Curated vernacular alias ("mouse" -> mmu, not "gray mouse lemur").
    alias = _ALIASES.get(key)
    if alias:
        hit = df.filter(pl.col("_code_lc") == alias)
        if not hit.is_empty():
            return _row_to_info(hit.row(0, named=True))

    # 1. KEGG organism code — the common case, checked first.
    hit = df.filter(pl.col("_code_lc") == key)
    if not hit.is_empty():
        return _row_to_info(hit.row(0, named=True))

    # 2. T-number / taxonomy id.
    hit = df.filter(
        (pl.col("ktax_id").str.to_lowercase() == key) | (pl.col("tax_id") == raw)
    )
    if not hit.is_empty():
        return _row_to_info(hit.row(0, named=True))

    # 3. Exact scientific or common name.
    hit = df.filter((pl.col("_sci_lc") == key) | (pl.col("_com_lc") == key))
    if not hit.is_empty():
        return _row_to_info(hit.row(0, named=True))

    # 4. "Homo sapiens (human)" — strip the parenthesised alias and retry both.
    m = re.match(r"^(.*?)\s*\(([^)]*)\)\s*$", raw)
    if m:
        for part in (m.group(1), m.group(2)):
            part = part.strip().lower()
            if not part:
                continue
            hit = df.filter((pl.col("_sci_lc") == part) | (pl.col("_com_lc") == part))
            if not hit.is_empty():
                return _row_to_info(hit.row(0, named=True))

    # 5. Genus abbreviation: "E. coli", "M. musculus".
    m = re.match(r"^([a-z])\.?\s+(\w+)", key)
    if m:
        rx = f"^{m.group(1)}[a-z]+ {re.escape(m.group(2))}"
        hit = df.filter(pl.col("_sci_lc").str.contains(rx))
        if not hit.is_empty():
            return _row_to_info(hit.row(0, named=True))

    # 6. Unique substring match.
    hit = df.filter(
        pl.col("_sci_lc").str.contains(re.escape(key), literal=False)
        | pl.col("_com_lc").str.contains(re.escape(key), literal=False)
    )
    if hit.height == 1:
        return _row_to_info(hit.row(0, named=True))
    if hit.height > 1:
        # Rank rather than taking an arbitrary first row: prefer a
        # word-boundary hit, then the shortest name (the least-qualified
        # organism is the one a bare vernacular name refers to), then the
        # lowest T-number for a stable tiebreak.
        word = hit.filter(
            pl.col("_com_lc").str.contains(rf"\b{re.escape(key)}\b")
            | pl.col("_sci_lc").str.contains(rf"\b{re.escape(key)}\b")
        )
        pool = word if not word.is_empty() else hit
        ranked = pool.with_columns(
            (pl.col("_sci_lc").str.len_chars() + pl.col("_com_lc").str.len_chars())
            .alias("_len")
        ).sort(["_len", "ktax_id"])
        return _row_to_info(ranked.row(0, named=True))

    raise SpeciesNotFoundError(raw, search_organisms(raw, n=5)["kegg_code"].to_list())


#: Backwards-compatible alias for the v2.x name.
def kegg_species_code(species: str = "hsa") -> str:
    """
    KEGG organism code for *species*, matching R's ``kegg.species.code()``.

    R returns the bare code string, so this does too.  Use
    :func:`get_species_code` when the full :class:`SpeciesInfo` record is
    wanted.

    >>> kegg_species_code("human")
    'hsa'
    """
    return get_species_code(species).kegg_code


# ---------------------------------------------------------------------------
# Optional refresh from KEGG
# ---------------------------------------------------------------------------

def refresh_organism_table(force: bool = False) -> Path:
    """
    Re-download KEGG's organism list and store it as a user-level override.

    Purely optional: the bundled table already covers every organism KEGG
    published at build time.  This exists so a long-lived install can pick up
    newly sequenced organisms without a package upgrade.

    Returns the path to the refreshed table.

    Raises NetworkError when KEGG is unreachable — the bundled table stays in
    use, so lookups keep working.
    """
    from .cache import http_get

    text = http_get(f"{KEGG_BASE}/list/organism", ttl=0 if force else 86400)

    rows: list[dict] = []
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        ktax, code, name = parts[0].strip(), parts[1].strip(), parts[2].strip()
        m = re.match(r"^(.*?)\s*\(([^)]*)\)\s*$", name)
        sci, com = (m.group(1), m.group(2)) if m else (name, "")
        rows.append({
            "ktax_id": ktax, "tax_id": "", "kegg_code": code,
            "scientific_name": sci, "common_name": com,
            "entrez_gnodes": "1", "kegg_geneid": "", "ncbi_geneid": "",
            "ncbi_proteinid": "", "uniprot": "",
        })

    if not rows:
        raise NetworkError("KEGG organism list came back empty; keeping bundled table.")

    # Preserve the richer cross-reference columns of the bundled table.
    bundled = {r["kegg_code"]: r for r in _load_table().to_dicts()}
    for r in rows:
        old = bundled.get(r["kegg_code"])
        if old:
            for c in ("tax_id", "entrez_gnodes", "kegg_geneid", "ncbi_geneid",
                      "ncbi_proteinid", "uniprot"):
                r[c] = old.get(c, "") or ""

    out = _user_table_path()
    cols = list(rows[0])
    body = "\n".join(["\t".join(cols)] + ["\t".join(str(r[c]) for c in cols) for r in rows])
    with gzip.open(out, "wt", encoding="utf-8") as fh:
        fh.write(body)

    _load_table.cache_clear()
    get_species_code.cache_clear()
    return out


def organism_count() -> int:
    """Number of organisms in the active table."""
    return _load_table().height


@lru_cache(maxsize=1)
def bioconductor_orgdb_table() -> pl.DataFrame:
    """
    The ``bods`` table: KEGG code <-> Bioconductor OrgDb package <-> default
    gene id type.  Used to pick a sensible default ID type per species.
    """
    with gzip.open(_BODS, "rt", encoding="utf-8") as fh:
        return pl.read_csv(fh.read().encode(), separator="\t", infer_schema_length=0)


def default_gene_idtype(species_code: str) -> str:
    """Return the id type KEGG uses for *species_code* ("ENTREZ" or "KEGG")."""
    tbl = bioconductor_orgdb_table()
    hit = tbl.filter(pl.col("kegg_code") == species_code)
    if not hit.is_empty():
        raw = str(hit.row(0, named=True)["id_type"]).upper()
        # bods stores Bioconductor's short codes ("eg", "tair", "orf");
        # translate to the gene_idtype names pathview() accepts.
        return {"EG": "ENTREZ"}.get(raw, raw)
    try:
        return "ENTREZ" if get_species_code(species_code).entrez_gnodes else "KEGG"
    except SpeciesNotFoundError:
        return "KEGG"

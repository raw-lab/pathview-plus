"""
bundled.py
Loading the reference tables that ship inside the wheel.

One reader, used everywhere.  Four modules previously carried their own copy
of this three-line helper, and each copy had to independently get the parsing
options right; one of them did not, because compound names such as
``5'-deoxy-5'-(methylthio)adenosine`` contain quote characters that a CSV
reader will happily treat as field delimiters.  These files are strict TSV
with no quoting, so quoting is disabled explicitly.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl

DATA_DIR = Path(__file__).parent / "data"


def read_bundled_tsv(name: str | Path) -> pl.DataFrame:
    """
    Read a bundled gzipped TSV as all-text columns.

    Parameters
    ----------
    name:
        File name inside ``pathview/data`` or an absolute path.
    """
    path = Path(name)
    if not path.is_absolute():
        path = DATA_DIR / path
    if not path.exists():
        raise FileNotFoundError(
            f"Bundled data file missing: {path.name}. The package was built "
            "without its data directory; reinstall pathview-plus."
        )
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return pl.read_csv(
            fh.read().encode(),
            separator="\t",
            quote_char=None,        # these files are unquoted TSV
            infer_schema_length=0,  # identifiers stay text: see the CLI note
        )


def bundled_files() -> dict[str, int]:
    """Bundled data files and their sizes in bytes — useful for diagnostics."""
    return {p.name: p.stat().st_size for p in sorted(DATA_DIR.glob("*.tsv.gz"))}

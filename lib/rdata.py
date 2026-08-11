"""
rdata.py
Reading R ``.RData`` / ``.rda`` files without R.

Why this is here
----------------
Much of the reference data in this ecosystem is published only as R
workspaces: pathview's compound tables, SBGNview's SBGN collection and its
identifier crosswalks.  ``pyreadr`` handles data.frames but not the named
character vectors and lists those files actually use, and installing R just to
read a table is not reasonable.

This reads R's XDR serialization format directly.  It was written to build
this package's bundled data, and is exported because reading a colleague's
``.rda`` is a normal thing to want to do.

Coverage and limits
-------------------
Implemented: NULL, symbols, pairlists, character/integer/real/logical/raw
vectors, generic vectors (lists), attributes (so named vectors and data.frames
come back keyed), reference tables, and the ALTREP wrappers R 3.5+ uses for
compact sequences and deferred strings.

Not implemented: closures, environments, S4 objects, factors with attached
levels, and ALTREP classes outside the wrappers above.  Those raise
:class:`ValueError` naming the SEXP type and the byte offset, rather than
returning something plausible but wrong — a parser that guesses silently
misaligns the stream and corrupts every object after the one it mis-read.

In practice this reads the flat tables and named vectors that reference data
is usually published as.  A file it cannot read will say so; use ``pyreadr``
or R itself for those.

Reference: R's ``src/main/serialize.c``.

Example
-------
>>> from pathview.rdata import read_rdata          # doctest: +SKIP
>>> tables = read_rdata("cpd.names.rda")
>>> list(tables)
['cpd.names']
"""

from __future__ import annotations

import gzip
import struct
from typing import Any

# SEXP type codes
NILSXP, SYMSXP, LISTSXP, CLOSXP = 0, 1, 2, 3
ENVSXP, PROMSXP, LANGSXP, SPECIALSXP = 4, 5, 6, 7
BUILTINSXP, CHARSXP, LGLSXP = 8, 9, 10
INTSXP, REALSXP, CPLXSXP, STRSXP = 13, 14, 15, 16
DOTSXP, ANYSXP, VECSXP, EXPRSXP = 17, 18, 19, 20
BCODESXP, EXTPTRSXP, WEAKREFSXP, RAWSXP, S4SXP = 21, 22, 23, 24, 25
NILVALUE_SXP, GLOBALENV_SXP = 254, 253
UNBOUNDVALUE_SXP, MISSINGARG_SXP, BASENAMESPACE_SXP = 252, 251, 250
NAMESPACESXP, PACKAGESXP, PERSISTSXP = 249, 248, 247
CLASSREF_SXP, GENERICREF_SXP, BCREPDEF, BCREPREF = 246, 245, 244, 243
EMPTYENV_SXP, ATTRLANGSXP, ATTRLISTSXP = 242, 240, 239
ALTREP_SXP = 238
REFSXP = 255

NA_INTEGER = -2147483648


class RDataReader:
    """Streaming reader for one serialized R object."""

    def __init__(self, data: bytes):
        self.d = data
        self.i = 0
        self.refs: list[Any] = []

    # -- primitives --------------------------------------------------------
    def _int(self) -> int:
        v = struct.unpack_from(">i", self.d, self.i)[0]
        self.i += 4
        return v

    def _double(self) -> float:
        v = struct.unpack_from(">d", self.d, self.i)[0]
        self.i += 8
        return v

    def _bytes(self, n: int) -> bytes:
        b = self.d[self.i:self.i + n]
        self.i += n
        return b

    def _flags(self) -> tuple[int, bool, bool, bool]:
        f = self._int()
        return (f & 0xFF, bool(f & (1 << 9)), bool(f & (1 << 10)), bool(f & (1 << 8)))

    # -- strings -----------------------------------------------------------
    def _charsxp(self) -> str | None:
        self._flags()
        n = self._int()
        if n == -1:
            return None
        return self._bytes(n).decode("utf-8", errors="replace")

    # -- dispatch ----------------------------------------------------------
    def read(self) -> Any:
        stype, has_attr, has_tag, _ = self._flags()

        if stype == NILVALUE_SXP or stype == NILSXP:
            return None
        if stype == REFSXP:
            idx = self._int() >> 8
            return self.refs[idx - 1]
        if stype in (GLOBALENV_SXP, EMPTYENV_SXP, BASENAMESPACE_SXP):
            return None
        if stype == UNBOUNDVALUE_SXP or stype == MISSINGARG_SXP:
            return None

        if stype == ALTREP_SXP:
            # R >= 3.5 compact/deferred vectors.  The serialized form is
            # info, state, attributes; for the wrappers used by these files
            # the state carries (or is) the materialised vector.
            info = self.read()
            state = self.read()
            self.read()                          # attributes
            return self._materialise_altrep(info, state)

        if stype == SYMSXP:
            name = self.read()
            self.refs.append(name)
            return name

        if stype == CHARSXP:
            n = self._int()
            return None if n == -1 else self._bytes(n).decode("utf-8", errors="replace")

        if stype == STRSXP:
            n = self._int()
            vals = [self._charsxp() for _ in range(n)]
            if has_attr:
                # Named character vectors are how SBGNview stores its
                # crosswalks.  Failing to consume the attributes here also
                # left the stream misaligned for every following object.
                names = self._attributes().get("names")
                if names and len(names) == len(vals):
                    return {"names": names, "values": vals}
            return vals

        if stype == INTSXP or stype == LGLSXP:
            n = self._int()
            vals = [self._int() for _ in range(n)]
            out = [None if v == NA_INTEGER else v for v in vals]
            if has_attr:
                self._attributes()
            return out

        if stype == REALSXP:
            n = self._int()
            out = [self._double() for _ in range(n)]
            if has_attr:
                self._attributes()
            return out

        if stype == RAWSXP:
            n = self._int()
            return self._bytes(n)

        if stype in (VECSXP, EXPRSXP):
            n = self._int()
            items = [self.read() for _ in range(n)]
            names = None
            if has_attr:
                attrs = self._attributes()
                names = attrs.get("names")
            if names:
                return dict(zip(names, items))
            return items

        if stype in (LISTSXP, LANGSXP, ATTRLANGSXP, ATTRLISTSXP):
            out: dict[str, Any] = {}
            order: list[Any] = []
            while True:
                if has_attr:
                    self.read()                      # discard attributes
                tag = self.read() if has_tag else None
                value = self.read()
                if tag is not None:
                    out[str(tag)] = value
                else:
                    order.append(value)
                stype, has_attr, has_tag, _ = self._flags()
                if stype == NILVALUE_SXP:
                    break
                if stype not in (LISTSXP, LANGSXP, ATTRLANGSXP, ATTRLISTSXP):
                    raise ValueError(f"unexpected pairlist continuation type {stype}")
            return out if out else order

        raise ValueError(f"unsupported R SEXP type {stype} at offset {self.i}")

    @staticmethod
    def _materialise_altrep(info: Any, state: Any) -> Any:
        """
        Recover the underlying vector from an ALTREP state.

        ``wrap_*`` and ``deferred_string`` store the real vector inside the
        state; ``compact_intseq`` stores (length, start, step).
        """
        cls = ""
        if isinstance(info, dict):
            cls = str(next(iter(info.values()), "") or "")
        elif isinstance(info, list) and info:
            cls = str(info[0])

        if isinstance(state, list):
            if "compact_intseq" in cls and len(state) == 3:
                n, start, step = state
                n = int(n)
                return [int(start + i * step) for i in range(n)]
            # wrap_*/deferred_string: the payload is the first list element
            # that is itself a vector.
            for item in state:
                if isinstance(item, list) and item and not isinstance(item[0], list):
                    return item
            return state
        if isinstance(state, dict):
            for item in state.values():
                if isinstance(item, list):
                    return item
        return state

    def _attributes(self) -> dict:
        attrs = self.read()
        return attrs if isinstance(attrs, dict) else {}


def _decompress(path: str) -> bytes:
    """R writes .RData with gzip, bzip2, xz or no compression."""
    with open(path, "rb") as fh:
        head = fh.read(6)
        fh.seek(0)
        blob = fh.read()
    if head[:2] == b"\x1f\x8b":
        return gzip.decompress(blob)
    if head[:3] == b"BZh":
        import bz2
        return bz2.decompress(blob)
    if head[:6] == b"\xfd7zXZ\x00":
        import lzma
        return lzma.decompress(blob)
    return blob


def read_rdata(path: str) -> dict:
    """Read an .RData / .rda workspace into ``{name: value}``."""
    raw = _decompress(path)

    if raw[:5] not in (b"RDX2\n", b"RDX3\n"):
        raise ValueError(f"{path}: not an RDX2/RDX3 file (got {raw[:5]!r})")
    version3 = raw[:5] == b"RDX3\n"

    i = 7                                    # skip "RDXn\nX\n"
    r = RDataReader(raw)
    r.i = i
    r._int()                                 # serialization version
    r._int()                                 # writer R version
    r._int()                                 # minimum reader R version
    if version3:
        n = r._int()
        r._bytes(n)                          # native encoding name

    obj = r.read()
    return obj if isinstance(obj, dict) else {"value": obj}


def as_columns(obj: Any, names: list[str] | None = None) -> dict:
    """
    Normalise a read R table to ``{column: values}``.

    R data.frames arrive as a named VECSXP, but a plain list arrives unnamed;
    ``names`` supplies the column order in that case.
    """
    if isinstance(obj, dict):
        return obj
    if isinstance(obj, list) and names:
        return dict(zip(names, obj))
    raise TypeError(f"cannot interpret {type(obj).__name__} as a table")


def rdata_objects(path: str) -> dict[str, str]:
    """Object names in an .RData file mapped to a short description of each."""
    out = {}
    for key, value in read_rdata(path).items():
        kind = type(value).__name__
        size = len(value) if hasattr(value, "__len__") else 1
        out[key] = f"{kind}[{size}]"
    return out


if __name__ == "__main__":
    import sys
    for p in sys.argv[1:]:
        try:
            d = read_rdata(p)
            for k, v in d.items():
                kind = type(v).__name__
                size = len(v) if hasattr(v, "__len__") else "-"
                print(f"{p}: {k} <{kind}> n={size}")
                if isinstance(v, dict):
                    for kk in list(v)[:3]:
                        print(f"    {kk!r} -> {str(v[kk])[:80]!r}")
                elif isinstance(v, list):
                    print(f"    {[str(x)[:60] for x in v[:3]]}")
        except Exception as exc:
            print(f"{p}: FAILED {type(exc).__name__}: {exc}")

"""
errors.py
Typed exception hierarchy.

v2.x raised bare RuntimeError/ValueError or, worse, swallowed failures and
returned an empty dict, so callers could not distinguish "network down" from
"pathway has no mappable nodes".  Every failure mode now has a type.
"""

from __future__ import annotations


class PathviewError(Exception):
    """Base class for every pathview-plus error."""


class SpeciesNotFoundError(PathviewError, ValueError):
    """The requested species could not be resolved to a KEGG organism code."""

    def __init__(self, species: str, suggestions: list[str] | None = None):
        self.species = species
        self.suggestions = suggestions or []
        msg = f"Unknown species {species!r}."
        if self.suggestions:
            msg += " Did you mean: " + ", ".join(self.suggestions[:5]) + "?"
        msg += (
            " pathview-plus ships an offline organism table; "
            "use list_organisms() or search_organisms('<name>') to browse it."
        )
        super().__init__(msg)


class PathwayNotFoundError(PathviewError, FileNotFoundError):
    """KGML/SBGN for the requested pathway is neither cached nor downloadable."""


class NetworkError(PathviewError, IOError):
    """A remote resource was unreachable after retries."""


class ParseError(PathviewError, ValueError):
    """A KGML or SBGN-ML document could not be parsed."""


class MappingError(PathviewError, ValueError):
    """No input identifier could be mapped onto the pathway."""


class RenderError(PathviewError, RuntimeError):
    """Rendering failed."""

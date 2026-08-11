# Contributing to Pathview-plus

Thanks for considering a contribution.

## Ground rules

**No stubs.** A function that warns and returns `None` is worse than no
function: it advertises a capability that does not exist, and callers cannot
tell the difference. If something is not implemented, say so in
`lib/parity.py` and raise a specific error naming the alternative.

`tests/test_parity.py` enforces this. It fails the build if a feature is marked
`full` without resolving to a real, importable API, and it walks the exported
surface looking for warn-only functions.

**Every bug fix gets a regression test** that fails against the version that
had the bug. See `tests/test_bug_regressions.py` — each test names the bug, its
symptom and its fix, so a reviewer can trace claim → test → fix.

**The test suite makes no network calls.** CI runs with `PATHVIEW_OFFLINE=1`.
If your change needs data, add a fixture under `tests/fixtures/` or bundle a
table under `lib/data/`. Tests that genuinely need the network go behind
the `PATHVIEW_TEST_NETWORK=1` guard, as in `tests/test_sbgn_collection.py`.

## Setup

```bash
git clone https://github.com/raw-lab/pathview-plus
cd pathview-plus
pip install -e ".[dev,docs]"

pytest -q
ruff check lib tests
cd docs && sphinx-build -b html . _build/html
```

## Before opening a pull request

1. `pytest -q` passes — currently 316 tests
2. `ruff check lib tests` is clean
3. New behaviour has a test; fixed bugs have a regression test
4. `lib/parity.py` updated if the change adds or removes a capability
5. `CHANGELOG.md` updated
6. Docstrings explain *why*, not just *what* — especially where the code
   departs from the obvious approach

## Building the bundled data

`lib/data/` is generated from upstream sources, not hand-edited. The
organism table comes from KEGG, the compound tables from R pathview, and the
SBGN index and crosswalks from the SBGNview / SBGNhub project via
`lib/rdata.py`. If you regenerate them, note the source and date in the
commit message.

## Areas that would help most

1. **ENSEMBL crosswalks** — two SBGNview files use an ALTREP variant
   `lib/rdata.py` does not decode, so ENSEMBL has no offline route to
   SBGN glyph ids
2. **Edge routing** — A* pathfinding around obstacles for splines
3. **SBGN glyph shapes** — a fuller SBGN-PD shape vocabulary
4. **Layout** — automatic relayout for maps with no coordinates
5. **Performance** — parallel batch rendering

Please open an issue before starting anything large, so we can agree on the
approach first.

## Reporting bugs

A good report includes: the pathway id and species, the render mode, a minimal
data frame that reproduces it, and the full traceback. `pathview-plus info`
output helps too.

## Licence

Contributions are accepted under the project's licence, CC BY-NC 4.0.

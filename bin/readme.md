# scripts here

## `pathview-cli.py`

Direct-run launcher, kept for continuity. It works straight from a clone
without installing anything:

```bash
python bin/pathview-cli.py --version
python bin/pathview-cli.py render 00020 --gene-data expr.csv
```

After `pip install`, the same interface is available as two commands on your
PATH, and those are the preferred way in:

```bash
pathview-plus render 00020 --gene-data expr.csv
pathview-cli   render 00020 --gene-data expr.csv   # alias
```

The implementation lives in `lib/cli.py` (`pathview.cli`); this file only
resolves the import path and calls `main()`.

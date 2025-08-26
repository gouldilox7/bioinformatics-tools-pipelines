# Repo Sanitizer — normalize file references to `study.<ext>`

Install:
```bash
python -m pip install --upgrade pyyaml
```

Dry-run (no changes):
```bash
python study.py --root /path/to/project --config patterns.yml --dry-run
```

Apply (content changes in-place + optional notebook stripping):
```bash
python study.py --root /path/to/project --config patterns.yml --apply --strip-notebooks
```

### What this does
- Finds file references ending with these extensions: study.csv, .tsv, .txt, .job, .sh, .batch, .R, .py`
- Rewrites the **basename** (just the last path component, before the extension) to `study`, preserving the directory and extension.
  - Example: `../data/zapata.sample.v1.csv` → `../data/study.csv`
  - Example: `"study.sh"` → `"study.sh"`
- **Content-only**: It does not rename files on disk (to avoid ambiguity/collisions). If you want on-disk renames too, we can enable a safe mode with a one-to-one map.

### Notes
- Works with multi-dot names (`foo.bar.csv`) and both `/` and `\` path separators.
- Skips changing names that are already `study.<ext>`.
- You can extend or shrink the `extensions` list in `patterns.yml`.

# Run several screens

Install Fit first. Save a JSON list such as `screens.json` beside your input
files; file paths are resolved relative to that JSON file:

```json
[
  {"name": "screen-a", "count_table": "a/counts.tsv", "design_matrix": "a/design.tsv"},
  {"name": "screen-b", "count_table": "b/counts.tsv", "design_matrix": "b/design.tsv"}
]
```

From the repository root:

```bash
python packages/crisprworks-fit/examples/batch.py \
  --screens screens.json --out-dir results --dry-run
python packages/crisprworks-fit/examples/batch.py \
  --screens screens.json --out-dir results --threads 4 --permutation-round 10
```

The runner checks all names, input paths and output directories before starting
any screen. It runs screens sequentially, stops on failure and writes each
screen's ordinary Fit outputs and provenance manifest in its own directory.
Existing screen directories require an explicit `--overwrite`. Commands are
executed as argument lists, without shell interpretation. Only the core MLE
options shown here are set; normalization, dispersion and efficiency behavior
retain the CLI defaults. Use the CLI directly for more elaborate designs or
efficiency updates.

Ten permutation rounds follow the upstream suggestion and cost more than the
two-round default used in our benchmark records. Choose settings appropriate
to your analysis rather than treating benchmark settings as scientific advice.

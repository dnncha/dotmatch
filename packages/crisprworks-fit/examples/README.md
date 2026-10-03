# CRISPRWorks Fit examples

## Complete joint analysis without assembling input files

From an installation of Fit:

```bash
crisprworks-fit demo --model joint --out-dir fit-joint-demo/
```

This offline synthetic example generates counts, a `Sample/Condition/Control`
map and 32 predeclared neutral genes. It writes joint effects, gene-control
calibration, named conditions, and `demo.json` with input/output hashes.
It requires a new directory. The default BY adjustment may make no calls
because finite control counts limit resolution. A successful demo checks the
software workflow, not biological accuracy.

Follow the [Count-to-Fit tutorial](https://dotmatch.readthedocs.io/en/latest/tutorials/crispr-fit-first-run.html)
for your own counts, model selection and status interpretation.

## Run several MLE screens

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
options shown here are set. Screens without a control-gene file retain the
CLI normalization defaults; dispersion retains the defaults. Guide-efficiency
updates are opt-in through `--update-efficiency`. Use the CLI directly for more
elaborate designs.

To use independently chosen nonessential genes for normalization and the
permutation null, add a `"control_gene": "a/nonessential-controls.txt"` field
to a screen. The file must contain one gene name per line, without a header.
The runner preflights its path and sets `--norm-method control`. Select finite
tails for the batch with `--permutation-pvalues finite`; the default `legacy`
mode retains upstream compatibility. The ordinary per-screen manifests record
the control file's hash, selected inference mode and null resolution.

To retain named, full-precision effects for later whole-gene calibration:

```bash
python packages/crisprworks-fit/examples/batch.py \
  --screens screens.json --out-dir calibrated-inputs \
  --permutation-pvalues finite --update-efficiency --write-fit-details
crisprworks-fit calibrate \
  --fit-details calibrated-inputs/screen-a/screen.fit-details.json \
  --control-gene a/nonessential-controls.txt -n calibrated-results/screen-a
```

The last command assumes the manifest's `screen-a` uses that independently
chosen control file. Full-precision outputs record MLE coefficient labels and
effect units. Existing outputs are preserved by choosing new result roots.

Ten permutation rounds follow the upstream suggestion and cost more than the
two-round default used in our benchmark records. Choose settings appropriate
to your analysis rather than treating benchmark settings as scientific advice.

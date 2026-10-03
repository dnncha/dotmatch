# From CRISPR counts to gene effects

CRISPRWorks Fit analyzes guide-count tables after counting QC. Install it
separately from CRISPRWorks Count. This tutorial covers a complete offline
example, the Count-to-Fit handoff and the choices that affect interpretation.
Fit is an experimental source package; it has no published PyPI release.

## 1. Install and run a complete example

From the DotMatch repository root, using Python 3.10 or newer:

```bash
python3 -m venv .venv-fit
. .venv-fit/bin/activate
python3 -m pip install ./packages/crisprworks-fit
crisprworks-fit --version
crisprworks-fit demo --model joint --out-dir fit-joint-demo/
```

Installation needs a C++ compiler and `make` for MAGeCK2. A C compiler builds
Fit's optional native MLE loop. The joint model uses NumPy and does not require
the native loop. The package's [pinned container and dependency instructions](https://github.com/dnncha/dotmatch/tree/main/packages/crisprworks-fit)
provide another installation route.

The joint demo runs locally without downloading data. It generates 48 genes
with four guides each, two biological-condition groups and a shared replicated
baseline. The generator declares 32 neutral control genes before fitting. It
uses whole-gene calibration with the default global BY adjustment. These are
synthetic software checks, not validated biological hits.

| File in `fit-joint-demo/` | What to inspect |
| --- | --- |
| `counts.tsv` | Raw guide counts with `sgRNA`, `Gene` and six sample columns |
| `samples.tsv` | Sample, condition and matched-baseline declarations for joint inference |
| `controls.txt` | One declared synthetic control gene per line |
| `design.tsv` | A labeled design for trying the same counts with MLE |
| `screen.joint.tsv` | Gene/condition effects and descriptive posterior uncertainty |
| `screen.joint-details.json` | Full-precision effects, guide efficacies, sample mappings, versions and hashes |
| `screen.calibration.tsv` | Named conditions, empirical tails, adjusted values and testing status |
| `screen.calibration.json` | Null strata, control counts, p-value resolution and testing family |
| `demo.json` | Completion marker and hashes of the example's inputs and outputs |

The joint demo refuses an existing output directory. If a run fails, inspect
its diagnostics and use a new directory. `demo.json` appears only after the
fit and calibration finish. The original MLE demo remains available as
`crisprworks-fit demo --out-dir fit-mle-demo/`; `--threads`, `--backend` and
`--kernel` select MLE execution behavior.

## 2. Hand off reviewed Count outputs

First complete the [Count tutorial](crispr-count-first-run.md). Review sample
identity, assignment settings, guide representation, zero counts and replicate
QC before fitting. Count's small installation demo is too small for joint
variance smoothing; use the joint example above to check Fit installation.

A completed Count assay writes
`crispr-screen/assay_out/counts.mageck.tsv`. This has the format Fit expects:

```text
sgRNA	Gene	baseline1	baseline2	treated1	treated2
guide_a	GENE_A	120	115	24	29
guide_b	GENE_A	90	104	21	23
```

This two-guide snippet illustrates columns, not a runnable joint dataset.
Use raw counts and the recorded gene annotations. Preserve exact guide and
sample IDs, including capitalization. Do not supply log-transformed counts or
quietly combine technical lanes into biological replicates. Joint counts must
be finite and nonnegative; guide IDs must be unique and gene IDs nonempty.
Control-gene files contain gene IDs, not guide IDs. Missing controls are not
silently replaced by other genes.

## 3. Choose the effect model before examining hits

| Analysis | Choose | Effect scale |
| --- | --- | --- |
| Fit a labeled experimental design while retaining MAGeCK2 conventions | `mle` | Natural-log beta coefficients for non-baseline design columns |
| Learn relative guide efficacy across replicated conditions using the same library | `joint` | Log2 relative-abundance effects, with one shared relative efficacy per guide |
| Test frozen effects against a whole-gene control null | `calibrate` after either fit | Empirical p-values and adjusted values; effects are not refitted |

Choose controls from independent assay knowledge or a predeclared annotation
set. They must represent null genes after the entire normalization and fitting
procedure. A gene neutral in another cell line need not be neutral in your
treatment. Selecting controls because their observed effects are small can
invalidate the null.

### MLE with a design matrix

Create `design.tsv` with sample labels matching the count columns:

```text
Samples	baseline	treatment
baseline1	1	0
baseline2	1	0
treated1	1	1
treated2	1	1
```

The first design column is all ones. The first row is a baseline with zeros
in the remaining columns. The design must identify the effects you intend to
estimate; include a documented full-rank design for additional factors.

```bash
crisprworks-fit mle \
  -k crispr-screen/assay_out/counts.mageck.tsv -d design.tsv \
  --control-gene controls.txt --norm-method control \
  --permutation-pvalues finite --update-efficiency \
  --seed 42 --permutation-round 10 --write-fit-details -n results/mle
```

The command writes upstream-compatible gene and guide summaries, a
`mle.crisprworks.json` provenance manifest and `mle.fit-details.json`.
Full-precision details preserve coefficient names from the design; an unlabeled
inline design uses `beta_1`, `beta_2`, and so on. MLE's default permutation
mode is `legacy`; this example explicitly chooses inclusive add-one tails.
More permutations improve resolution, not the validity of a pooled-guide null.
The round count is an example, not an assay-independent recommendation.

### Joint inference with a sample map

Create a tab-separated `samples.tsv` with the exact header below:

```text
Sample	Condition	Control
baseline1	BASE	BASE
baseline2	BASE	BASE
treated1	TREATED	BASE
treated2	TREATED	BASE
```

Map every count-table sample exactly once. Baselines map to themselves;
treatment conditions name a baseline condition and require at least two
replicates. A baseline may have one replicate. Add further replicated
conditions to share guide-efficacy estimates, and use separate baseline groups
when the experiment requires them. The complete table needs at least 64 guides.

```bash
crisprworks-fit joint \
  -k crispr-screen/assay_out/counts.mageck.tsv --sample-map samples.tsv \
  --control-gene controls.txt -n results/joint
```

The outputs are `joint.joint.tsv` and `joint.joint-details.json`. Conditions
are sorted by name and recorded explicitly. Negative effects indicate relative
depletion; positive effects indicate relative enrichment. Relative guide
efficacies are continuous model weights, not editing probabilities. Posterior
standard deviations and z scores describe the fitted model; they are not
calibrated biological hit probabilities.

The model adapts JACKS's default guide-efficacy inference. It does not provide
copy-number correction, pretrained guide references or an automatic test of
the difference between two fitted conditions. MLE beta and joint log2 effects
should be interpreted on their recorded scales rather than compared as equal
numbers.

## 4. Calibrate and read testing status

Calibrate the full-precision JSON, rather than a rounded summary table:

```bash
crisprworks-fit calibrate \
  --fit-details results/joint.joint-details.json \
  --control-gene controls.txt --adjust by --family global -n results/calibrated
```

For MLE, replace `--fit-details` with `results/mle.fit-details.json`. The same
procedure can be rerun on the offline example using
`fit-joint-demo/screen.joint-details.json` and `fit-joint-demo/controls.txt`.

The calibrator compares whole control-gene effects within exact guide-count
strata. It defaults to at least 20 controls per stratum. The MLE permutation
limit excludes MLE genes with 40 or more guides; successful joint fits are not
excluded by this MLE-specific limit. With B controls, the smallest directional empirical
p-value is `1/(B+1)`, before multiple-testing adjustment. A supported stratum
does not guarantee enough resolution to discover hits.

| `status` | Interpretation |
| --- | --- |
| `tested` | A supported finite score was calibrated; this does not itself declare a hit |
| `training_control` | Used for null calibration; excluded from the testing family and assigned p=q=1 |
| `insufficient_controls` | Too few controls in this guide-count stratum; p=q=1 |
| `invalid_fit` | Nonfinite tested effect; p=q=1 and an empty score field |
| `invalid_control_fit` | A nonfinite control invalidates this condition's null stratum; p=q=1 |
| `skipped_fit` | MLE record excluded by the configured permutation guide-count limit; p=q=1 |

`condition` is a one-based coefficient index; `condition_label` preserves the
design or joint-map label when the input contains names. Older detail files may
contain indices alone, so retain their original design or sample map.
`p_negative`/`q_negative` describe the lower tail and
`p_positive`/`q_positive` the upper tail. Use `p_two`/`q_two` when choosing
direction after inspecting effects. Each alternative has its own declared
testing family.

BY is the default adjustment for arbitrary dependence of valid marginal
p-values. `--adjust bh` needs independence or suitable positive dependence.
The default global family includes every non-control gene/condition hypothesis;
`--family condition` declares a separate family per condition. Choose these
settings before examining calls. Neither BH nor BY repairs unrepresentative
controls or invalid fitted null scores.

No calls can be a legitimate result of limited null resolution and conservative
adjustment. Do not change controls or adjustment just to obtain more hits.
Check the manifest's control counts, family size and statuses alongside the
biological controls and analysis plan.

## 5. Retain provenance and the scope of the evidence

Keep Count's assignment settings and QC with the raw-count table, design or
sample map, control declarations, software versions and complete Fit outputs.
The manifests record hashes for verification. Fit commands write outputs at
their supplied prefixes; choose new prefixes to preserve prior results. Only
the joint demo automatically refuses an existing output directory.

The [scientific validation report](https://github.com/dnncha/dotmatch/blob/main/packages/crisprworks-fit/benchmarks/SCIENTIFIC-VALIDATION.md)
separates numerical parity, runtime, essential-gene ranking and simulated error
calibration. Current joint comparisons use curated control-enriched subsets
from three cell lines and two libraries. They do not establish overall SOTA,
whole-genome FDR or independent biological accuracy. Count QC, synthetic demos
and a faster MLE fit do not establish those claims either.

Joint and calibration tables are complete only when their companion JSON
manifest exists and its output hash matches the table. Both files are staged
before publication; interrupted publication can leave a table without a
manifest. Overlapping joint/calibration publishers at the same prefix are
rejected before replacement; choose distinct prefixes to retain both runs.
Hidden `.publish.lock` files remain in place and are not completion markers;
the OS releases their locks when a process exits. Do not delete a lock file
while publishers are running. This is not a guarantee of durability across
filesystem or machine failure.

Joint details now include each gene's `fit_diagnostics`: iteration count,
termination reason, convergence flag and final absolute stopping-statistic
change. `joint --max-iterations 50 --tolerance 0.1` shows the defaults; both
options require positive values and tolerance must be finite. A convergence
flag means the inherited stopping rule was met, not that posterior accuracy
was independently established. Calibration rejects explicit nonconvergence;
legacy files without diagnostics remain readable but have no convergence check.

When fit details record a companion summary, calibration requires that summary
to exist at its recorded path and match its SHA-256 hash. Keep both files when
moving or archiving a joint result and update the recorded path when relocating
it. Legacy standalone detail files without a companion record remain supported.
The hash check verifies summary consistency; it does not authenticate the
fit-detail contents or validate their biological assumptions.

Calibration manifests now include a `resolution` diagnostic. It applies the
same BH/BY adjustment and global/condition family to optimistic p-value floors,
while retaining unsupported, skipped and invalid hypotheses at p=1. It reports
minimum attainable directional/two-sided q-values and the number of hypotheses
whose bounds reach the diagnostic threshold. `--fdr-alpha` selects that threshold
(default 0.05, valid range `(0, 1]`); it changes neither p/q values nor the family.
The CLI warns when no directional or two-sided discovery is attainable. No calls
in that situation cannot establish absence of effects. An optimistic bound does
not predict power or establish exchangeability of the controls.

Joint and calibration parse their inputs from immutable in-memory bytes and
record hashes of exactly those bytes. They recheck the input files and output
aliases immediately before publishing. A changed or missing input aborts
publication, preserving an earlier completed bundle. These checks and advisory
locks do not prevent unrelated software from subsequently altering files.

MLE also fits from private input snapshots and stages the entire required output
bundle. Failed fitting or staging leaves previous completed outputs intact. It
uses the same publication lock and manifest-last policy; interruption during
replacement can leave incomplete files without a completion manifest. A
successful rerun without `--write-fit-details` removes earlier optional details.
Do not treat an older complete manifest as proof that a failed new attempt
succeeded: check the command exit status and verify the recorded output hashes.

The MLE count preflight rejects inconsistent rows, empty/duplicate columns,
empty guide/gene identifiers, duplicate guide IDs and nonfinite or negative
counts before fitting. CSV snapshots preserve their `.csv` suffix for upstream
delimiter selection; parsing retains upstream literal TSV/CSV behavior.

# Package integrity audit — 6 September 2026

Baseline: released 0.5.0, main commit
`11d159fa1648365f2a4e96917b483c33aa5d9fe7`, tree
`8cdb16395656d31655a08e0449ef29720bec7679`.

## Why this pass focuses on integrity

An auditable count table is only useful when installation, parsing, label
identity and downstream conversion preserve its meaning. This audit followed
those boundaries rather than adding another top-level product surface. The
native distance kernel and default assignment semantics are not changed.

## Reproduced findings and changes

| Boundary | Baseline behaviour | Change |
| --- | --- | --- |
| Python thresholds | `k=1.9` rounds; `2**32` wraps through ctypes | Integer-only, explicit representable range, no booleans |
| Experimental posterior | NaN/Infinity priors propagate or cause incidental errors; equal candidates can be called unique at a low cutoff | Finite priors, log-space normalization and explicit ambiguous ties |
| Raw count import | Missing values become zero; 9007199254740993 rounds through float | Strict widths and lossless integer parsing; missing counts require explicit zero |
| Dataframe input | Documented read DataFrames fail with AttributeError; column/label selection is inconsistent | Named and explicit original-axis column selection, exact label lengths, no null stringification |
| AnnData count bridge | Permissive fallback and numeric guessing bypass malformed-table errors or sample selection | One strict count reader; sparse int64 matrix; source sample names/order and metadata retained |
| Cell-feature bridge | Dense pivot; cell labels inferred from read IDs; zero-count cells dropped | Sparse CSR, explicit cells, optional fixed axes, per-cell observation accounting |
| Matcher lifetime | `close()` can free an index while a ctypes call uses it | Per-instance reentrant lock; deterministic concurrent-close regression |
| Native runtime lookup | CWD/arbitrary ancestor lookup and silent fallback after a bad override | Packaged/source-layout discovery and authoritative explicit overrides |
| Version identity | Unrelated ancestor pyproject can override package version | Only the actual dotmatch source layout can supply the source version |
| Python FASTQ | Duplicated parsers disagree on malformed identifiers and quality | Shared strict four-line reader, contextual errors, original-byte digest semantics |
| Native target tables | Fixed 16 KiB rows, limited fields, no compressed/quoted CSV support | Bounded dynamic gzip-aware parser, aligned headers, malformed input rejected |
| Container | Native-only entrypoint bypasses sensitivity/assay/agent; compiler and source in runtime | Packaged unified dispatcher, portable packaging build, multi-stage runtime |

The baseline had many passing tests. Their existence did not cover these input
boundaries. The new tests reproduce the specific failures instead of relaxing
existing expected scientific results. The pre-existing feature-matrix test is
updated only to assert the deliberate sparse return type; its counts and labels
are unchanged.

## Acceptance evidence

Local environment: Python 3.13 on Linux with pandas/NumPy/SciPy. Polars and
AnnData were unavailable locally; their integration assertions must be executed
in the new Scientific Python CI jobs, not represented by local skips. Docker
build/runtime checks likewise execute in the dedicated Linux x86-64 and ARM64
jobs. Refer to the PR's exact-commit checks and attached test XML for final
counts and completion state.

The tests include: independent scalar Levenshtein dynamic programming at lengths
0/1/8/20/31/32/33/63/64/65/127/128/129; indexed and unindexed policy comparisons;
invalid numeric parameters; posterior ties; real native/Python CSV, gzip and BOM
agreement; long-row and corrupt-input failures; complete count conservation;
exact values above 2**53; sparse H5AD round trips; zero-count axes; and a
50,000-by-50,000 sparse matrix with one nonzero. The prior exhaustive eight-base
sensitivity oracle remains unchanged.

The container check uses the built image, network disabled, a read-only root
filesystem, a caller-owned output mount and an unprivileged UID. It exercises
both native and Python commands and checks that compiler/git tools and the
source checkout are absent from the runtime. Images in this audit workflow are
not pushed and the existing public 0.5.0 image is not overwritten.

## Compatibility and limits

See the Unreleased changelog and the dataframe guide for intentional rejection
of malformed formerly-coerced values. Native one-column libraries now use
Python's `target_0` convention. Named IDs remain unchanged. Library IDs/genes
starting with a literal quote are rejected to protect downstream TSV identity.
The native target parser accepts at most 128 columns and 1 MiB per logical
single-line row; it rejects rather than truncates an oversized row.

This is not a complete certification of every command. Native high-throughput
FASTQ parsing remains a separate implementation; fuzzing and cross-reader
validation of every native fast path are further work. Low-level APIs compare
literal bytes, and the experimental posterior model is not empirically
calibrated. Real external scientific accuracy and end-to-end competitor
throughput need independent experimental/benchmark evidence. No new such
claims are made here. Historical benchmark artifacts and the native scoring
kernel are unchanged.

## Next decision-relevant work

Prioritize native FASTQ/fuzz coverage and output crash-recovery before more
formats. Extend policy-sensitivity reports to representative public full-size
libraries with matched resource/semantic baselines. Validate one maintained
upstream workflow adapter with an independent lab. Introduce quality-based
assignment or cell-level inference only with their corresponding labelled
validation, not by relabelling geometry or per-read agreement as accuracy.

## Implementation references

- AnnData supports sparse CSR/CSC matrices and annotated axes:
  https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html
- ctypes native calls can run without the interpreter lock:
  https://docs.python.org/3/library/ctypes.html
- Docker build stages separate build tools from the runtime:
  https://docs.docker.com/build/building/multi-stage/
- GitHub-hosted runner platform labels:
  https://docs.github.com/en/actions/reference/runners/github-hosted-runners

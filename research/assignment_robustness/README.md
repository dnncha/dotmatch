# AR-001: assignment robustness in pooled CRISPR reads

**Read the [findings](RESULTS.md).** Five public-data investigations are complete, including independently checked counting over two full Yusa archives (20,394,663 reads). Both historical Yusa and Brunello discrepancies reproduce exactly. The principal finding is that guide counting can mix read records with matching-offset events, and the gene-annotation pattern of positional ambiguity differs markedly between the tested libraries.

This is a research branch, not a product release or a completed gene-discovery paper. No biological hit list has been generated. The 600,000-read pilot, historical prefixes and full Yusa archives overlap; their denominators must not be added as independent data.

## Study records

| File | Purpose |
|---|---|
| [PROTOCOL.md](PROTOCOL.md) | Original exploratory pilot protocol, prior observations disclosed |
| [AMENDMENTS.md](AMENDMENTS.md) | Failed-run audit trail and pre-execution full-Yusa extension |
| [BRUNELLO_FORENSICS.md](BRUNELLO_FORENSICS.md) | Prespecified historical Brunello reconstruction |
| [RESULTS.md](RESULTS.md) | Technical results, units, validation and limitations |
| `pilot.py` | Discovery/evaluation-separated six-sample controlled comparison |
| `forensics.py` | Exact reconstruction of the historical Yusa event counts |
| `full_yusa.py` | Complete-archive checksum, record and guide-count verification |
| `brunello_forensics.py` | Historical Brunello multi-offset event reconstruction |
| `library_aliases.py` | Post-hoc exact short-shift overlap graph and tests |
| `summarize.py` | Offline verification of immutable artifacts and descriptive summaries |
| `test_summary.py` | Corrupt/missing evidence, unsafe path and trace regression tests |

All native production code remains pinned to baseline `11d159fa1648365f2a4e96917b483c33aa5d9fe7`. Execution commits and run links are in RESULTS.md. The comparator is guide-counter 0.1.3, not an unpinned latest install.

## Reproduce the calculations from retained outputs

Download the five successful evidence ZIPs from the run links in RESULTS.md, or use the complete evidence bundle supplied with the study handoff. Keep each ZIP in a separate extraction directory; do not overwrite files from one run with another run's files.

With the directory names used below, no network, raw sequencing files or compiled DotMatch executable is needed to validate and regenerate the summary:

```bash
python research/assignment_robustness/summarize.py \
  --pilot pilot/study-output \
  --yusa-prefix yusa-prefix/historical-yusa \
  --brunello-prefix brunello-prefix/brunello-forensic-output \
  --overlap overlap/library-alias-output \
  --full-yusa full-yusa/full-yusa-output \
  --out regenerated-summary
```

The output directory must not exist. The summary verifies 190 manifest-listed files, rejects incomplete runs, checks reference hashes across lanes and recomputes the read/event, gene-annotation and overlap statistics from retained outputs. `results.json`, `prefix_accounting.tsv` and `fixed_window_pilot.tsv` are generated. This is an offline reconstruction of reported statistics, not a substitute for rerunning raw-read assignment.

## Reproduce the raw-data experiments

Use a clean checkout of this research branch or the recorded execution commits. Python 3.11, a supported C++ build environment and Cargo are required for this path. No paid infrastructure is needed by the scripts. They download public data and write local files; the full Yusa input archives are about 573 MB compressed. Each output directory must be new.

```bash
python -m pip install --no-deps .
cargo install guide-counter --version 0.1.3 --locked

python research/assignment_robustness/pilot.py test
python research/assignment_robustness/full_yusa.py test
python research/assignment_robustness/library_aliases.py test
python research/assignment_robustness/test_summary.py

python research/assignment_robustness/pilot.py run --out pilot-output
python research/assignment_robustness/forensics.py yusa --out yusa-prefix-output
python research/assignment_robustness/brunello_forensics.py --out brunello-prefix-output
python research/assignment_robustness/library_aliases.py --out overlap-output
python research/assignment_robustness/full_yusa.py --out full-yusa-output
```

The fixed-window pilot intentionally does not provide a complete Brunello extraction workflow. Its window selection is a controlled diagnostic. The historical multi-offset lane uses different documented extraction semantics and is analysed separately.

Runtime and absolute paths can differ between reruns, so complete manifests need not be byte-identical. Source data hashes, count vectors, denominators and defined biological annotations are the important invariants. A changed upstream reference or FASTQ must be treated as a changed input rather than silently accepted as replication of the original source.

## Completed versus outstanding

Completed: source inspection, six-sample oracle validation, both historical discrepancy reconstructions, full Yusa confirmation, both complete-library structural audits and offline evidence verification. Two initial research-harness failures are preserved and explained rather than erased.

Outstanding: full Brunello sample processing, assay-aware extraction suitable for gene analysis, a locked replicated downstream cohort and analysis plan, matched-depth sensitivity controls, held-out studies, independent biological validation, a broader novelty review, external scientific review and a publication deposit.

The next biological cohort candidate is the source study's modified-tracr Brunello A375 dropout experiment. The retained ENA metadata maps RepA to SRR8297837 + SRR8297836, RepB to SRR8297839 + SRR8297838, RepC to SRR8297841 + SRR8297840, and plasmid to SRR8297997. Those paired runs are technical lanes of their respective biological samples, not six biological replicates. The shared plasmid is not three independent controls. These are candidate inputs, not a newly approved full-screen analysis plan.

Before downstream hit calling, freeze extraction and sample metadata against the source study, pin the downstream software and settings, set inclusion/filtering and multiplicity rules, retain null results and separate pilot-informed studies from held-out validation. A new significance threshold crossing is not a biological discovery without adequate independent support.

## Evidence stewardship

Public Actions artifacts are retained for 30 days and are scheduled to expire on 6 October 2026. The handoff bundle preserves the successful outputs, source and failed-run evidence. No raw FASTQs are redistributed in git or the evidence artifacts. Source terms still apply to the input reads and library files; this project does not relicense them. A durable archival deposit with a DOI remains required before manuscript submission.

The report is a technical Markdown artifact; native DotMatch sensitivity HTML reports are retained inside the experiment outputs. Exact tables are used for audit lookup, rather than presenting these two selected libraries as an adequately sampled population. This research was implemented with AI assistance and has not yet received independent human scientific review. Final authorship and disclosure should be approved before submission. Do not cite this work as peer reviewed.

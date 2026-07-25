# Lab Evaluation and Handoff

This guide is for a bioinformatics team or core facility evaluating DotMatch on
a known-target sequencing assay. It is a local technical evaluation protocol,
not a clinical validation protocol and not a substitute for assay-specific
controls.

## Before processing study data

Confirm all of the following with the assay owner:

- the target library is the intended revision and contains only the expected
  fixed-window sequences;
- the FASTQ read and orientation are known, along with the start and length of
  the target window;
- the permitted edit radius and ambiguity policy have been chosen deliberately;
- sample identifiers and FASTQ paths have been checked against the run sheet;
- the intended downstream consumer of the output is known. For CRISPR counts,
  DotMatch writes a MAGeCK-style matrix but does not perform screen statistics.

Create a reviewable project from the release package:

```bash
dotmatch assay new crispr \
  --library guides.csv \
  --reads-dir fastqs/ \
  --out crispr_evaluation/

cd crispr_evaluation
dotmatch assay check assay.toml
```

`assay new` samples the input reads to propose an extraction window. Review
`inference_report.json`, the target-library audit, and the generated
`assay.toml`. Keep `status = "draft"` until a qualified reviewer has confirmed
the configuration, then set it to `ready`.

## Run and review

```bash
./run.sh
dotmatch assay handoff assay.toml
```

`./run.sh` uses `dotmatch assay start`: it runs preflight, assignment, target
audit, QC, and validation. Open these files in order:

1. `assay_out/reliability_report.html`
2. `assay_out/sample_qc.tsv` and `assay_out/crispr_qc.html` for CRISPR runs
3. `assay_out/assay_report.html`
4. `assay_out/counts.mageck.tsv` or the primary count/demultiplexing output
5. `assay_out/methods.md` and `assay_out/CITATION.bib`

Do not treat a `passed` DotMatch reliability verdict as proof of biological
validity. It means the configured target safety and software QC rules passed.
Review positive/negative controls, sample identity, sequencing-run metrics,
replicate agreement, and downstream analysis according to the local assay
protocol.

## Handoff package

`dotmatch assay handoff assay.toml` writes `assay_out/handoff/` without copying
raw reads. The bundle is suitable for an internal technical review or a
workflow-maintainer evaluation:

- `README_FOR_REVIEW.md` states the review order and boundary;
- `handoff_manifest.json` records the configuration, verdict, input file sizes,
  and input/output SHA-256 hashes;
- `SHA256SUMS` verifies the copied review files;
- `review/` contains reports, QC tables, primary outputs, methods, citation,
  and software-version records.

In the controlled workspace containing the FASTQs, recompute each input hash
in `handoff_manifest.json` before approving the handoff. Do not send raw reads
or identifiers outside the approved data-handling route merely to make a
DotMatch review bundle.

## Decision record

For each evaluated assay, record:

| Item | Record |
| --- | --- |
| Assay and target-library revision | Name, source, checksum, and owner |
| Read extraction | Read, orientation, start, length, and rationale |
| Assignment rules | Metric, edit radius, ambiguity policy, and handling of ambiguous reads |
| Input identity | Sample-sheet revision and FASTQ checksums |
| QC outcome | Reliability verdict, findings reviewed, controls, and any exceptions |
| Output recipient | Count matrix/report location and downstream analysis owner |
| Software record | DotMatch version, native version, `methods.md`, and `CITATION.bib` |

Keep this record with the project or laboratory notebook. It makes a later
rerun auditable without claiming that DotMatch replaces the remainder of the
assay or analysis workflow.

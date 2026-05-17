import assert from "node:assert/strict";
import test from "node:test";

import {
  buildCountAssayToml,
  buildDemuxAssayToml,
  buildPairAssayToml,
  candidateSummary,
  manifestSummaryColumns,
  parseManifestSummary,
  promoteDraftSpec,
  qcStatus
} from "../src/lib/assayModel.js";

test("buildCountAssayToml creates a ready AssaySpec with workspace-relative paths", () => {
  const toml = buildCountAssayToml({
    status: "ready",
    assayType: "crispr",
    targets: "inputs/guides.tsv",
    samples: [
      { id: "control", fastq: "reads/control.fastq.gz" },
      { id: "treated", fastq: "reads/treated.fastq.gz" }
    ],
    outDir: "runs/crispr",
    start: 23,
    length: 19,
    k: 1,
    metric: "levenshtein"
  });

  assert.match(toml, /schema_version = 1/);
  assert.match(toml, /status = "ready"/);
  assert.match(toml, /mode = "count"/);
  assert.match(toml, /assay_type = "crispr"/);
  assert.match(toml, /targets = "inputs\/guides.tsv"/);
  assert.match(toml, /id = "control"/);
  assert.match(toml, /fastq = "reads\/treated.fastq.gz"/);
  assert.match(toml, /\[extract\]\nstart = 23\nlength = 19/);
  assert.match(toml, /\[assignment\]\nk = 1\nmetric = "levenshtein"/);
});

test("buildDemuxAssayToml creates a demux AssaySpec", () => {
  const toml = buildDemuxAssayToml({
    status: "ready",
    assayType: "inline_barcode",
    barcodes: "inputs/barcodes.tsv",
    reads: "reads/pooled.fastq.gz",
    outDir: "runs/demux",
    start: 0,
    length: 8,
    k: 1,
    metric: "hamming"
  });

  assert.match(toml, /mode = "demux"/);
  assert.match(toml, /barcodes = "inputs\/barcodes.tsv"/);
  assert.match(toml, /reads = "reads\/pooled.fastq.gz"/);
  assert.match(toml, /\[extract\]\nstart = 0\nlength = 8/);
});

test("buildPairAssayToml creates a pair-count AssaySpec", () => {
  const toml = buildPairAssayToml({
    status: "ready",
    assayType: "generic",
    leftTargets: "inputs/left.tsv",
    rightTargets: "inputs/right.tsv",
    reads: "reads/pairs.fastq.gz",
    outDir: "runs/pair",
    leftStart: 0,
    leftLength: 8,
    rightStart: 20,
    rightLength: 8,
    k: 1,
    metric: "hamming"
  });

  assert.match(toml, /mode = "pair-count"/);
  assert.match(toml, /left_targets = "inputs\/left.tsv"/);
  assert.match(toml, /right_targets = "inputs\/right.tsv"/);
  assert.match(toml, /\[left\]\nstart = 0\nlength = 8/);
  assert.match(toml, /\[right\]\nstart = 20\nlength = 8/);
});

test("buildCountAssayToml escapes TOML strings instead of injecting raw values", () => {
  const toml = buildCountAssayToml({
    status: "draft",
    assayType: "generic",
    targets: "inputs/guide\"bad.tsv",
    samples: [{ id: "bad\nsample", fastq: "reads/a.fastq" }],
    outDir: "runs/out",
    start: 0,
    length: 8,
    k: 0,
    metric: "hamming"
  });

  assert.match(toml, /status = "draft"/);
  assert.match(toml, /targets = "inputs\/guide\\"bad.tsv"/);
  assert.match(toml, /id = "bad\\nsample"/);
  assert.doesNotMatch(toml, /id = "bad\nsample"/);
});

test("promoteDraftSpec changes only the top-level draft status", () => {
  const promoted = promoteDraftSpec('schema_version = 1\nstatus = "draft"\nmode = "count"\n');

  assert.equal(promoted, 'schema_version = 1\nstatus = "ready"\nmode = "count"\n');
});

test("parseManifestSummary keeps documented columns deterministic", () => {
  const tsv =
    "schema_version\tmode\tassay_type\tstatus\tnative_version\tautopsy_triggered\twarning_count\tproduction_warning_count\tsample_count\tprimary_report\tmanifest\n" +
    "1\tcount\tcrispr\tready\tdotmatch 0.1.0\tfalse\t0\t1\t2\tassay_report.html\tassay_manifest.json\n";

  const rows = parseManifestSummary(tsv);

  assert.deepEqual(manifestSummaryColumns, [
    "schema_version",
    "mode",
    "assay_type",
    "status",
    "native_version",
    "autopsy_triggered",
    "warning_count",
    "production_warning_count",
    "sample_count",
    "primary_report",
    "manifest"
  ]);
  assert.equal(rows.length, 1);
  assert.equal(rows[0].mode, "count");
  assert.equal(rows[0].production_warning_count, "1");
});

test("qcStatus flags samples that cross automatic autopsy thresholds", () => {
  assert.equal(qcStatus({ assignment_rate: 0.92, ambiguous_rate: 0.01, no_match_rate: 0.04, invalid_rate: 0 }), "pass");
  assert.equal(qcStatus({ assignment_rate: 0.79, ambiguous_rate: 0.01, no_match_rate: 0.04, invalid_rate: 0 }), "review");
  assert.equal(qcStatus({ assignment_rate: 0.9, ambiguous_rate: 0.08, no_match_rate: 0.04, invalid_rate: 0 }), "review");
  assert.equal(qcStatus({ assignment_rate: 0.9, ambiguous_rate: 0.01, no_match_rate: 0.16, invalid_rate: 0 }), "review");
  assert.equal(qcStatus({ assignment_rate: 0.9, ambiguous_rate: 0.01, no_match_rate: 0.04, invalid_rate: 0.03 }), "review");
});

test("candidateSummary ranks candidate windows by confidence and separation", () => {
  const summary = candidateSummary({
    chosen: { start: 23, length: 19, assignment_rate: 0.91, ambiguous_rate: 0.01, no_match_rate: 0.08, score: 0.9, score_margin: 0.29 },
    candidates: [
      { start: 21, length: 19, assignment_rate: 0.7, ambiguous_rate: 0.45, no_match_rate: 0.1, score: 0.25 },
      { start: 22, length: 19, assignment_rate: 0.62, ambiguous_rate: 0.01, no_match_rate: 0.37, score: 0.61 },
      { start: 23, length: 19, assignment_rate: 0.91, ambiguous_rate: 0.01, no_match_rate: 0.08, score: 0.9, score_margin: 0.29 }
    ],
    warnings: []
  });

  assert.equal(summary.confidence, "high");
  assert.equal(summary.separation, 0.29);
  assert.equal(summary.rejected[0].start, 22);
});

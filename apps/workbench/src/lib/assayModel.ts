export type AssayStatus = "ready" | "draft";
export type AssayMetric = "hamming" | "levenshtein";
export type QcReviewStatus = "pass" | "review";

export type CountSampleInput = {
  id: string;
  fastq: string;
};

export type CountAssayInput = {
  status: AssayStatus;
  assayType: "crispr" | "feature_barcode" | "inline_barcode" | "amplicon_panel" | "oligo_adapter" | "generic";
  targets: string;
  samples: CountSampleInput[];
  outDir: string;
  start: number;
  length: number;
  k: number;
  metric: AssayMetric;
};

export type DemuxAssayInput = {
  status: AssayStatus;
  assayType: CountAssayInput["assayType"];
  barcodes: string;
  reads: string;
  outDir: string;
  start: number;
  length: number;
  k: number;
  metric: AssayMetric;
};

export type PairAssayInput = {
  status: AssayStatus;
  assayType: CountAssayInput["assayType"];
  leftTargets: string;
  rightTargets: string;
  reads: string;
  outDir: string;
  leftStart: number;
  leftLength: number;
  rightStart: number;
  rightLength: number;
  k: number;
  metric: AssayMetric;
};

export type QcRates = {
  assignment_rate: number;
  ambiguous_rate: number;
  no_match_rate: number;
  invalid_rate: number;
};

export type InferenceCandidate = {
  start: number;
  length: number;
  assignment_rate: number;
  ambiguous_rate: number;
  no_match_rate: number;
  score?: number;
  score_margin?: number;
};

export type InferenceReport = {
  chosen: InferenceCandidate;
  candidates: InferenceCandidate[];
  warnings: string[];
};

export type CandidateReview = {
  confidence: "high" | "medium" | "low";
  separation: number;
  rejected: InferenceCandidate[];
};

export const manifestSummaryColumns = [
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
] as const;

export type ManifestSummaryColumn = (typeof manifestSummaryColumns)[number];
export type ManifestSummaryRow = Record<ManifestSummaryColumn, string>;

export function buildCountAssayToml(input: CountAssayInput): string {
  const lines = [
    "schema_version = 1",
    `status = ${tomlString(input.status)}`,
    'mode = "count"',
    `assay_type = ${tomlString(input.assayType)}`,
    `targets = ${tomlString(input.targets)}`,
    "",
    ...input.samples.flatMap((sample) => [
      "[[samples]]",
      `id = ${tomlString(sample.id)}`,
      `fastq = ${tomlString(sample.fastq)}`,
      ""
    ]),
    "[run]",
    `out_dir = ${tomlString(input.outDir)}`,
    "threads = 1",
    "",
    "[extract]",
    `start = ${input.start}`,
    `length = ${input.length}`,
    "",
    "[assignment]",
    `k = ${input.k}`,
    `metric = ${tomlString(input.metric)}`,
    'ambiguity_policy = "best"',
    'ambiguous = "discard"',
    "",
    "[outputs]",
    'format = "mageck"',
    "assignments = true",
    "ambiguous = true",
    "unmatched = true",
    ""
  ];

  return lines.join("\n");
}

export function buildDemuxAssayToml(input: DemuxAssayInput): string {
  return [
    "schema_version = 1",
    `status = ${tomlString(input.status)}`,
    'mode = "demux"',
    `assay_type = ${tomlString(input.assayType)}`,
    `barcodes = ${tomlString(input.barcodes)}`,
    `reads = ${tomlString(input.reads)}`,
    "",
    "[run]",
    `out_dir = ${tomlString(input.outDir)}`,
    "threads = 1",
    "",
    "[extract]",
    `start = ${input.start}`,
    `length = ${input.length}`,
    "",
    "[assignment]",
    `k = ${input.k}`,
    `metric = ${tomlString(input.metric)}`,
    'ambiguity_policy = "best"',
    'ambiguous = "discard"',
    "",
    "[outputs]",
    "assignments = true",
    "ambiguous = true",
    "unmatched = true",
    ""
  ].join("\n");
}

export function buildPairAssayToml(input: PairAssayInput): string {
  return [
    "schema_version = 1",
    `status = ${tomlString(input.status)}`,
    'mode = "pair-count"',
    `assay_type = ${tomlString(input.assayType)}`,
    `left_targets = ${tomlString(input.leftTargets)}`,
    `right_targets = ${tomlString(input.rightTargets)}`,
    `reads = ${tomlString(input.reads)}`,
    "",
    "[run]",
    `out_dir = ${tomlString(input.outDir)}`,
    "threads = 1",
    "",
    "[left]",
    `start = ${input.leftStart}`,
    `length = ${input.leftLength}`,
    "",
    "[right]",
    `start = ${input.rightStart}`,
    `length = ${input.rightLength}`,
    "",
    "[assignment]",
    `k = ${input.k}`,
    `metric = ${tomlString(input.metric)}`,
    'ambiguity_policy = "best"',
    'ambiguous = "discard"',
    "",
    "[outputs]",
    "assignments = true",
    ""
  ].join("\n");
}

export function promoteDraftSpec(toml: string): string {
  if (/^status\s*=\s*"draft"\s*$/m.test(toml)) {
    return toml.replace(/^status\s*=\s*"draft"\s*$/m, 'status = "ready"');
  }
  return `status = "ready"\n${toml}`;
}

export function parseManifestSummary(tsv: string): ManifestSummaryRow[] {
  const rows = parseTsv(tsv);
  return rows.map((row) => {
    const normalized = Object.fromEntries(manifestSummaryColumns.map((column) => [column, row[column] ?? ""]));
    return normalized as ManifestSummaryRow;
  });
}

export function qcStatus(rates: QcRates): QcReviewStatus {
  if (rates.assignment_rate < 0.8) return "review";
  if (rates.ambiguous_rate > 0.05) return "review";
  if (rates.no_match_rate > 0.15) return "review";
  if (rates.invalid_rate > 0.02) return "review";
  return "pass";
}

export function candidateSummary(report: InferenceReport): CandidateReview {
  const chosenKey = candidateKey(report.chosen);
  const rejected = report.candidates
    .filter((candidate) => candidateKey(candidate) !== chosenKey)
    .sort((a, b) => candidateScore(b) - candidateScore(a));
  const secondBest = rejected[0] ? candidateScore(rejected[0]) : 0;
  const separation = roundRate(report.chosen.score_margin ?? candidateScore(report.chosen) - secondBest);
  const hasWarnings = report.warnings.length > 0;
  let confidence: CandidateReview["confidence"] = "low";

  if (!hasWarnings && report.chosen.assignment_rate >= 0.85 && separation >= 0.15) {
    confidence = "high";
  } else if (report.chosen.assignment_rate >= 0.7 && separation >= 0.05) {
    confidence = "medium";
  }

  return { confidence, separation, rejected };
}

export function parseTsv(text: string): Record<string, string>[] {
  const lines = text.trim().split(/\r?\n/).filter(Boolean);
  if (lines.length === 0) return [];
  const headers = lines[0].split("\t");
  return lines.slice(1).map((line) => {
    const values = line.split("\t");
    return Object.fromEntries(headers.map((header, index) => [header, values[index] ?? ""]));
  });
}

function tomlString(value: string): string {
  return `"${value.replace(/\\/g, "\\\\").replace(/"/g, '\\"').replace(/\n/g, "\\n").replace(/\r/g, "\\r")}"`;
}

function candidateKey(candidate: InferenceCandidate): string {
  return `${candidate.start}:${candidate.length}`;
}

function candidateScore(candidate: InferenceCandidate): number {
  return candidate.score ?? candidate.assignment_rate - candidate.ambiguous_rate;
}

function roundRate(value: number): number {
  return Math.round(value * 100) / 100;
}

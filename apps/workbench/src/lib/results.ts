import { parseTsv } from "./assayModel.js";

export type SampleQcRow = {
  sample_id: string;
  fastq: string;
  total_reads: number;
  assignment_rate: number;
  ambiguous_rate: number;
  no_match_rate: number;
  invalid_rate: number;
  [key: string]: string | number;
};

export type ManifestLike = {
  artifacts?: Record<string, string>;
  autopsy_artifacts?: Record<string, string>;
};

export type ArtifactGroup = {
  title: string;
  items: Array<{ label: string; path: string }>;
};

const numericQcFields = new Set([
  "total_reads",
  "valid_extracted_reads",
  "assigned_reads",
  "exact_reads",
  "k1_rescued_reads",
  "k1_sub_reads",
  "k1_ins_reads",
  "k1_del_reads",
  "ambiguous_reads",
  "no_match_reads",
  "invalid_reads",
  "assignment_rate",
  "exact_rate",
  "rescue_rate",
  "ambiguous_rate",
  "no_match_rate",
  "targets_observed",
  "zero_count_targets",
  "gini_index",
  "top_1pct_read_fraction",
  "candidates_verified"
]);

export function parseSampleQc(tsv: string): SampleQcRow[] {
  return parseTsv(tsv).map((row) => {
    const converted: Record<string, string | number> = {};
    for (const [key, value] of Object.entries(row)) {
      converted[key] = numericQcFields.has(key) ? Number(value || 0) : value;
    }
    return {
      sample_id: String(converted.sample_id ?? ""),
      fastq: String(converted.fastq ?? ""),
      total_reads: Number(converted.total_reads ?? 0),
      assignment_rate: Number(converted.assignment_rate ?? 0),
      ambiguous_rate: Number(converted.ambiguous_rate ?? 0),
      no_match_rate: Number(converted.no_match_rate ?? 0),
      invalid_rate: Number(converted.invalid_rate ?? 0),
      ...converted
    };
  });
}

export function safeDisplayText(value: string): string {
  return value
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#39;");
}

export function artifactGroups(manifest: ManifestLike): ArtifactGroup[] {
  const artifacts = manifest.artifacts ?? {};
  const autopsy = manifest.autopsy_artifacts ?? {};

  return [
    {
      title: "Primary report",
      items: itemList([["Assay report", artifacts.assay_report]])
    },
    {
      title: "Run provenance",
      items: itemList([
        ["Manifest JSON", artifacts.manifest],
        ["Manifest summary TSV", artifacts.manifest_summary],
        ["Normalized spec", artifacts.normalized_spec]
      ])
    },
    {
      title: "Quality control",
      items: itemList([
        ["Sample QC", artifacts.sample_qc],
        ["CRISPR QC report", artifacts.crispr_qc_report],
        ["CRISPR QC JSON", artifacts.crispr_qc],
        ["CRISPR QC summary", artifacts.crispr_qc_summary],
        ["Summary JSON", artifacts.summary],
        ["Native report", artifacts.report]
      ])
    },
    {
      title: "Audit and autopsy",
      items: itemList([
        ["Audit directory", artifacts.audit],
        ["Left audit", artifacts.left_audit],
        ["Right audit", artifacts.right_audit],
        ["Autopsy summary", autopsy.summary],
        ["Autopsy findings", autopsy.findings],
        ["Autopsy directory", autopsy.autopsy]
      ])
    },
    {
      title: "Workflow exports",
      items: itemList([
        ["Counts", artifacts.counts],
        ["Long counts", artifacts.target_counts_long],
        ["Pair counts", artifacts.pair_counts],
        ["Demuxed reads", artifacts.demuxed]
      ])
    }
  ];
}

function itemList(entries: Array<[string, string | undefined]>): ArtifactGroup["items"] {
  return entries
    .filter((entry): entry is [string, string] => Boolean(entry[1]))
    .map(([label, path]) => ({ label, path }));
}

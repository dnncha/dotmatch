import assert from "node:assert/strict";
import test from "node:test";

import { artifactGroups, parseSampleQc, safeDisplayText } from "../src/lib/results.js";

test("parseSampleQc converts numeric QC fields and preserves sample identifiers", () => {
  const rows = parseSampleQc(
    "sample_id\tfastq\ttotal_reads\tassignment_rate\tambiguous_rate\tno_match_rate\tinvalid_rate\n" +
      "sample_a\treads/a.fastq\t100\t0.94\t0.01\t0.05\t0\n"
  );

  assert.equal(rows.length, 1);
  assert.equal(rows[0].sample_id, "sample_a");
  assert.equal(rows[0].total_reads, 100);
  assert.equal(rows[0].assignment_rate, 0.94);
});

test("safeDisplayText escapes user-controlled file and warning values", () => {
  assert.equal(safeDisplayText('<img src=x onerror="alert(1)">'), "&lt;img src=x onerror=&quot;alert(1)&quot;&gt;");
});

test("artifactGroups exposes primary report, manifest, QC, audit, autopsy, and workflow exports", () => {
  const groups = artifactGroups({
    artifacts: {
      assay_report: "runs/out/assay_report.html",
      manifest: "runs/out/assay_manifest.json",
      manifest_summary: "runs/out/assay_manifest.summary.tsv",
      audit: "runs/out/audit",
      sample_qc: "runs/out/sample_qc.tsv",
      crispr_qc: "runs/out/crispr_qc.json",
      crispr_qc_report: "runs/out/crispr_qc.html",
      crispr_qc_summary: "runs/out/crispr_qc.summary.tsv"
    },
    autopsy_artifacts: {
      summary: "runs/out/autopsy/autopsy_summary.json",
      findings: "runs/out/autopsy/findings.tsv"
    }
  });

  assert.deepEqual(
    groups.map((group) => group.title),
    ["Primary report", "Run provenance", "Quality control", "Audit and autopsy", "Workflow exports"]
  );
  assert.equal(groups[0].items[0].path, "runs/out/assay_report.html");
  assert.equal(groups[2].items.some((item) => item.path.includes("crispr_qc.html")), true);
  assert.equal(groups[3].items.some((item) => item.path.includes("findings.tsv")), true);
});

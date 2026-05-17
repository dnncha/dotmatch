import { useMemo, useState } from "react";

import {
  buildCountAssayToml,
  buildDemuxAssayToml,
  buildPairAssayToml,
  candidateSummary,
  parseManifestSummary,
  qcStatus,
  type CountAssayInput
} from "./lib/assayModel.js";
import { artifactGroups, parseSampleQc, safeDisplayText } from "./lib/results.js";
import { doctor, readTextArtifact, runWorkbenchCommand, writeTextArtifact, type CommandResult, type DotmatchAction } from "./lib/workbenchApi.js";
import "./styles.css";

const defaultSpecPath = "assay.toml";

type JobRecord = CommandResult & { label: string };
type AssayMode = "count" | "demux" | "pair-count";

export default function App() {
  const [workspace, setWorkspace] = useState("");
  const [dotmatchPath, setDotmatchPath] = useState("");
  const [doctorStatus, setDoctorStatus] = useState<string>("Not checked");
  const [specPath, setSpecPath] = useState(defaultSpecPath);
  const [mode, setMode] = useState<AssayMode>("count");
  const [assayType, setAssayType] = useState<CountAssayInput["assayType"]>("crispr");
  const [targets, setTargets] = useState("examples/workflows/fixtures/crispr_library.csv");
  const [barcodes, setBarcodes] = useState("inputs/barcodes.tsv");
  const [leftTargets, setLeftTargets] = useState("inputs/left_targets.tsv");
  const [rightTargets, setRightTargets] = useState("inputs/right_targets.tsv");
  const [reads, setReads] = useState("reads/pooled.fastq.gz");
  const [sampleFastq, setSampleFastq] = useState("examples/workflows/fixtures/sample_a.fastq");
  const [outDir, setOutDir] = useState("workbench_out");
  const [start, setStart] = useState(0);
  const [length, setLength] = useState(8);
  const [k, setK] = useState(1);
  const [metric, setMetric] = useState<"hamming" | "levenshtein">("hamming");
  const [jobs, setJobs] = useState<JobRecord[]>([]);
  const [manifestText, setManifestText] = useState("");
  const [sampleQcText, setSampleQcText] = useState("");
  const [message, setMessage] = useState("");

  const toml = useMemo(
    () => {
      if (mode === "demux") {
        return buildDemuxAssayToml({
          status: "ready",
          assayType,
          barcodes,
          reads,
          outDir,
          start,
          length,
          k,
          metric
        });
      }
      if (mode === "pair-count") {
        return buildPairAssayToml({
          status: "ready",
          assayType,
          leftTargets,
          rightTargets,
          reads,
          outDir,
          leftStart: start,
          leftLength: length,
          rightStart: start + length,
          rightLength: length,
          k,
          metric
        });
      }
      return buildCountAssayToml({
        status: "ready",
        assayType,
        targets,
        samples: [{ id: "sample", fastq: sampleFastq }],
        outDir,
        start,
        length,
        k,
        metric
      });
    },
    [mode, assayType, barcodes, reads, outDir, start, length, k, metric, leftTargets, rightTargets, targets, sampleFastq]
  );

  const parsedManifest = useMemo(() => {
    try {
      return manifestText ? JSON.parse(manifestText) : {};
    } catch {
      return {};
    }
  }, [manifestText]);
  const qcRows = useMemo(() => (sampleQcText ? parseSampleQc(sampleQcText) : []), [sampleQcText]);
  const manifestSummaryRows = useMemo(() => {
    const summary = parsedManifest?.artifacts?.manifest_summary;
    if (!summary || !workspace) return [];
    return [];
  }, [parsedManifest, workspace]);
  void manifestSummaryRows;

  async function runDoctor() {
    setDoctorStatus("Checking DotMatch...");
    try {
      const report = await doctor(dotmatchPath || undefined);
      const failed = report.checks.filter((check) => check.exit_code !== 0);
      setDoctorStatus(
        failed.length === 0
          ? `Ready: ${report.dotmatch_path}`
          : `Found ${report.dotmatch_path}, but ${failed.length} doctor check failed`
      );
    } catch (error) {
      setDoctorStatus(String(error));
    }
  }

  async function saveSpec() {
    const path = await writeTextArtifact(workspace, specPath, toml);
    setMessage(`Wrote ${path}`);
  }

  async function runAction(label: string, action: DotmatchAction, args: string[]) {
    setMessage(`${label} running...`);
    const result = await runWorkbenchCommand({
      workspace,
      dotmatchPath: dotmatchPath || undefined,
      action,
      args,
      logName: label.toLowerCase().replace(/[^a-z0-9]+/g, "-")
    });
    setJobs((current) => [{ label, ...result }, ...current]);
    setMessage(`${label} finished with exit code ${result.exit_code}`);
  }

  function inferenceArgs(): string[] {
    const base = ["--mode", mode, "--assay-type", assayType, "--out", specPath, "--report", "inference_report.json", "--candidates", "inference_candidates.tsv"];
    if (mode === "demux") {
      return [...base, "--barcodes", barcodes, "--reads", reads];
    }
    if (mode === "pair-count") {
      return [...base, "--left-targets", leftTargets, "--right-targets", rightTargets, "--reads", reads];
    }
    return [...base, "--targets", targets, "--reads", sampleFastq];
  }

  async function refreshResults() {
    try {
      const manifest = await readTextArtifact(workspace, `${outDir}/assay_manifest.json`);
      setManifestText(manifest);
      const qc = await readTextArtifact(workspace, `${outDir}/sample_qc.tsv`);
      setSampleQcText(qc);
      const summaryPath = JSON.parse(manifest)?.artifacts?.manifest_summary;
      if (summaryPath) {
        const summary = await readTextArtifact(workspace, summaryPath);
        parseManifestSummary(summary);
      }
      setMessage("Loaded latest run artifacts");
    } catch (error) {
      setMessage(String(error));
    }
  }

  return (
    <main className="workbench-shell">
      <header className="topbar">
        <div>
          <p className="eyebrow">Local desktop app</p>
          <h1>DotMatch Workbench</h1>
        </div>
        <span className="local-badge">Data stays local</span>
      </header>

      <section className="panel setup-panel">
        <div>
          <h2>Workspace and DotMatch</h2>
          <p>Choose a local workspace. The runner only accepts relative paths inside that directory.</p>
        </div>
        <label>
          Workspace
          <input value={workspace} onChange={(event) => setWorkspace(event.target.value)} placeholder="/path/to/project" />
        </label>
        <label>
          DotMatch executable
          <input value={dotmatchPath} onChange={(event) => setDotmatchPath(event.target.value)} placeholder="PATH lookup" />
        </label>
        <button type="button" onClick={runDoctor}>Run doctor</button>
        <output>{doctorStatus}</output>
      </section>

      <section className="workbench-grid">
        <article className="panel">
          <h2>Design assay</h2>
          <div className="form-grid">
            <label>
              Mode
              <select value={mode} onChange={(event) => setMode(event.target.value as AssayMode)}>
                <option value="count">Count</option>
                <option value="demux">Demux</option>
                <option value="pair-count">Pair-count</option>
              </select>
            </label>
            <label>
              Assay type
              <select value={assayType} onChange={(event) => setAssayType(event.target.value as CountAssayInput["assayType"])}>
                <option value="crispr">CRISPR</option>
                <option value="feature_barcode">Feature barcode</option>
                <option value="inline_barcode">Inline barcode</option>
                <option value="amplicon_panel">Amplicon panel</option>
                <option value="oligo_adapter">Oligo adapter</option>
                <option value="generic">Generic</option>
              </select>
            </label>
            <label>Spec path<input value={specPath} onChange={(event) => setSpecPath(event.target.value)} /></label>
            <label>Targets<input value={targets} onChange={(event) => setTargets(event.target.value)} /></label>
            <label>Barcodes<input value={barcodes} onChange={(event) => setBarcodes(event.target.value)} /></label>
            <label>Left targets<input value={leftTargets} onChange={(event) => setLeftTargets(event.target.value)} /></label>
            <label>Right targets<input value={rightTargets} onChange={(event) => setRightTargets(event.target.value)} /></label>
            <label>Reads<input value={reads} onChange={(event) => setReads(event.target.value)} /></label>
            <label>Sample FASTQ<input value={sampleFastq} onChange={(event) => setSampleFastq(event.target.value)} /></label>
            <label>Output directory<input value={outDir} onChange={(event) => setOutDir(event.target.value)} /></label>
            <label>Start<input type="number" value={start} onChange={(event) => setStart(Number(event.target.value))} /></label>
            <label>Length<input type="number" value={length} onChange={(event) => setLength(Number(event.target.value))} /></label>
            <label>k<input type="number" value={k} onChange={(event) => setK(Number(event.target.value))} /></label>
            <label>
              Metric
              <select value={metric} onChange={(event) => setMetric(event.target.value as "hamming" | "levenshtein")}>
                <option value="hamming">Hamming</option>
                <option value="levenshtein">Levenshtein</option>
              </select>
            </label>
          </div>
          <div className="button-row">
            <button type="button" onClick={saveSpec} disabled={!workspace}>Save spec</button>
            <button type="button" onClick={() => runAction("Infer", "AssayInfer", inferenceArgs())} disabled={!workspace}>Infer</button>
            <button type="button" onClick={() => runAction("Check", "AssayCheck", [specPath])} disabled={!workspace}>Check</button>
            <button type="button" onClick={() => runAction("Plan", "AssayPlan", [specPath])} disabled={!workspace}>Plan</button>
            <button type="button" onClick={() => runAction("Run", "AssayRun", [specPath])} disabled={!workspace}>Run</button>
            <button type="button" onClick={() => runAction("Autopsy", "AssayAutopsy", [specPath, "--out-dir", `${outDir}/autopsy`])} disabled={!workspace}>Autopsy</button>
          </div>
          <pre className="code-preview"><code>{toml}</code></pre>
        </article>

        <article className="panel">
          <h2>Run status</h2>
          <p>{message || "No job has run in this session."}</p>
          <button type="button" onClick={refreshResults} disabled={!workspace}>Load results</button>
          <div className="job-list">
            {jobs.map((job) => (
              <details key={`${job.label}-${job.log_path}`}>
                <summary>{job.label}: exit {job.exit_code}</summary>
                <pre><code>{job.argv.join(" ")}</code></pre>
                <pre><code>{job.stdout || job.stderr}</code></pre>
              </details>
            ))}
          </div>
        </article>
      </section>

      <section className="workbench-grid">
        <article className="panel">
          <h2>QC review</h2>
          <div className="qc-grid">
            {qcRows.map((row) => (
              <div className={`qc-card ${qcStatus(row)}`} key={row.sample_id}>
                <strong>{row.sample_id}</strong>
                <span>{Math.round(row.assignment_rate * 100)}% assigned</span>
                <span>{Math.round(row.no_match_rate * 100)}% no match</span>
              </div>
            ))}
          </div>
        </article>

        <article className="panel">
          <h2>Artifacts</h2>
          <div className="artifact-list">
            {artifactGroups(parsedManifest).map((group) => (
              <section key={group.title}>
                <h3>{group.title}</h3>
                {group.items.length === 0 ? <p>No artifact recorded yet.</p> : null}
                {group.items.map((item) => (
                  <code key={`${group.title}-${item.path}`}>{safeDisplayText(`${item.label}: ${item.path}`)}</code>
                ))}
              </section>
            ))}
          </div>
        </article>
      </section>

      <section className="panel">
        <h2>Inference review</h2>
        <p>
          Candidate windows are shown after `assay infer` writes `inference_report.json`. Draft specs must be reviewed and marked ready before running.
        </p>
        <pre><code>{JSON.stringify(candidateSummary({
          chosen: { start, length, assignment_rate: 0, ambiguous_rate: 0, no_match_rate: 0 },
          candidates: [],
          warnings: ["Run inference to populate candidate windows."]
        }), null, 2)}</code></pre>
      </section>
    </main>
  );
}

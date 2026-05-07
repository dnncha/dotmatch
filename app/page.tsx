const useCases = [
  {
    title: "CRISPR guide counting",
    body: "Run count or crispr-count on FASTQ/FASTQ.gz and emit MAGeCK-compatible matrices, sample QC, target provenance, and summaries."
  },
  {
    title: "Inline barcode demultiplexing",
    body: "Split fixed-position single-end FASTQ/FASTQ.gz by barcode while keeping ambiguous and unmatched reads inspectable."
  },
  {
    title: "Classic BCL demultiplexing",
    body: "Process first-milestone classic per-cycle Illumina BCL run folders into sample FASTQ.gz and demultiplexing stats."
  },
  {
    title: "Audit, validation, and diagnosis",
    body: "Audit unsafe one-edit libraries, validate indexed runs against an exhaustive oracle, and inspect top unmatched reads."
  }
];

const proof = [
  ["Current scope", "known targets", "Short DNA barcodes, guides, primers, adapters, panels, and whitelist-style target sets. Not genome-scale alignment."],
  ["Assignment invariant", "indexed = exhaustive", "Indexed paths must match the native exhaustive assignment oracle for the same targets, k, and ambiguity policy."],
  ["Public CRISPR gate", "5x repeated", "Repeated 10k and 100k-record/sample MAGeCK/Yusa rows pass with count agreement and Edlib validation artifacts."],
  ["Open-core surface", "CLI + C + Python", "The C core, CLI, Python ctypes bindings, schemas, audit, validation, reports, and benchmarks are inspectable."]
];

const commands = [
  "dotmatch crispr-count --library guides.csv --samples samples.tsv --guide-start 23 --guide-length 19 --k 1 --metric levenshtein --indel-window 1 --out counts.mageck.tsv --summary qc.json",
  "dotmatch count --targets guides.csv --reads sample.fastq.gz --target-start 23 --target-length 19 --k 1 --metric levenshtein --indel-window 1 --report report.html --sample-qc sample_qc.tsv",
  "dotmatch demux --barcodes barcodes.tsv --reads pooled.fastq.gz --barcode-start 0 --barcode-length 8 --k 1 --metric hamming --out-dir demuxed --summary demux.qc.json",
  "dotmatch bcl-demux --run-folder 240101_RUN --sample-sheet SampleSheet.csv --out-dir bcl_demuxed --barcode-mismatches 1 --summary bcl.summary.json",
  "dotmatch audit --targets guides.tsv --k 1 --out-dir audit",
  "dotmatch inspect-unmatched --targets guides.tsv --reads sample.fastq.gz --target-start 23 --target-length 19 --k 1 --offset-window 2 --top 100 --out top_unmatched.tsv",
  "dotmatch validate --targets guides.tsv --reads sample.fastq.gz --target-start 23 --target-length 19 --k 1 --indel-window 1 --oracle edlib --sample 100000"
];

const throughputRows = [
  { label: "DotMatch exact k=0", value: 1143740, tone: "green" },
  { label: "DotMatch Hamming k=1", value: 331494, tone: "green" },
  { label: "guide-counter one mismatch", value: 194968, tone: "blue" },
  { label: "MAGeCK exact count", value: 92761, tone: "gray" },
  { label: "DotMatch Levenshtein k=1", value: 8836, tone: "green" }
] as const;

const memoryRows = [
  { label: "guide-counter one mismatch", value: 528.7, tone: "blue" },
  { label: "MAGeCK exact count", value: 158.9, tone: "gray" },
  { label: "DotMatch exact k=0", value: 28.7, tone: "green" },
  { label: "DotMatch Hamming k=1", value: 28.7, tone: "green" },
  { label: "DotMatch Levenshtein k=1", value: 27.5, tone: "green" }
] as const;

const candidateRows = [
  { label: "DotMatch Levenshtein verified/read", value: 2.822, tone: "green" },
  { label: "Exhaustive scan targets/read", value: 87437, tone: "blue" }
] as const;

const agreementRows = [
  { label: "DotMatch exact vs MAGeCK exact", value: 1.0, tone: "green" },
  { label: "DotMatch Hamming vs guide-counter", value: 0.942, tone: "blue" }
] as const;

const scalingRows = [
  { label: "guide-counter, 8 samples", value: 34800, tone: "blue" },
  { label: "guide-counter, 4 samples", value: 17400, tone: "blue" },
  { label: "guide-counter, 2 samples", value: 8700, tone: "blue" },
  { label: "DotMatch, 8 samples", value: 0, tone: "green" }
] as const;

export default function Home() {
  return (
    <main>
      <header className="site-header">
        <a className="brand" href="#top" aria-label="DotMatch home">
          <span className="brand-mark" />
          DotMatch
        </a>
        <nav aria-label="Primary navigation">
          <a href="#benchmarks">Benchmarks</a>
          <a href="#use-cases">Use cases</a>
          <a href="#comparison">Comparison</a>
          <a href="#workflow">Workflow</a>
        </nav>
        <a className="header-cta" href="#workflow">Get started</a>
      </header>

      <section id="top" className="hero">
        <div className="hero-copy">
          <h1>DotMatch</h1>
          <p className="hero-lede">
            Fast exact assignment for known short-DNA target sets.
          </p>
          <p className="hero-text">
            DotMatch is a C core, CLI, and Python binding for short barcode,
            guide, primer, adapter, and panel reads. It supports exact,
            Hamming, and k=1 Levenshtein assignment while keeping unique,
            ambiguous, no-match, and invalid outcomes explicit.
          </p>
          <div className="hero-actions">
            <a href="#benchmarks" className="button primary">
              Explore results
            </a>
            <a href="#workflow" className="button secondary">
              See commands
            </a>
          </div>
        </div>
        <div className="hero-panel" aria-label="DotMatch benchmark summary">
          <div className="panel-topline">
            <span>v0.1.0 scope</span>
            <span>known short-DNA targets</span>
          </div>
          <div className="metric-grid">
            <div>
              <strong>87,437</strong>
              <span>public MAGeCK/Yusa guides in the current claim-gate benchmark</span>
            </div>
            <div>
              <strong>0</strong>
              <span>Edlib validation mismatches across sampled public CRISPR reads</span>
            </div>
            <div>
              <strong>331k</strong>
              <span>reads/s DotMatch Hamming k=1 mean on repeated 100k-record/sample rows</span>
            </div>
            <div>
              <strong>28.7 MB</strong>
              <span>peak RSS for DotMatch Hamming k=1 repeated public CRISPR rows</span>
            </div>
          </div>
          <div className="sequence-rail" aria-hidden="true">
            {Array.from({ length: 64 }).map((_, i) => (
              <span key={i} className={i % 7 === 0 ? "hot" : i % 5 === 0 ? "cool" : ""} />
            ))}
          </div>
        </div>
      </section>

      <section id="benchmarks" className="section proof-section">
        <div className="section-heading">
          <h2>Current product and benchmark state.</h2>
          <p>
            The current project is an open-core short-DNA assignment engine, not
            a general aligner. The defensible public claim is known-target
            FASTQ assignment with explicit ambiguity semantics and reproducible
            benchmark artifacts.
          </p>
        </div>
        <div className="proof-grid">
          {proof.map(([label, value, detail]) => (
            <article className="proof-card" key={label}>
              <span>{label}</span>
              <strong>{value}</strong>
              <p>{detail}</p>
            </article>
          ))}
        </div>
        <div className="benchmark-grid">
          <article className="benchmark-card">
            <div className="chart-copy">
              <span className="card-label">Public CRISPR claim gate</span>
              <h3>Repeated Yusa rows are now the headline evidence.</h3>
              <p>
                Five repeated 100k-record/sample rows are checked in for
                DotMatch, MAGeCK, and guide-counter. The chart keeps exact,
                Hamming, and Levenshtein semantics in separate lanes.
              </p>
            </div>
            <HorizontalBarChart
              rows={throughputRows}
              unit="reads/s"
              axisLabel="Mean throughput, 100k records/sample, log scale"
              scale="log"
            />
          </article>

          <article className="benchmark-card">
            <div className="chart-copy">
              <span className="card-label">Candidate verification</span>
              <h3>Levenshtein k=1 verifies about three candidates per read.</h3>
              <p>
                Public Yusa k=1 Levenshtein rows average 2.822 verified
                candidates per read across 87,437 guides, then exact
                verification decides unique, ambiguous, or no-match status.
              </p>
            </div>
            <HorizontalBarChart
              rows={candidateRows}
              unit="checks/read"
              axisLabel="Work per read, log scale"
              scale="log"
            />
          </article>

          <article className="benchmark-card">
            <div className="chart-copy">
              <span className="card-label">Memory profile</span>
              <h3>CRISPR counting stays in tens of MB.</h3>
              <p>
                On the repeated public Yusa rows, DotMatch Hamming and exact
                count runs peak around 28.7 MB. guide-counter remains around
                528.7 MB on the same fixture.
              </p>
            </div>
            <HorizontalBarChart
              rows={memoryRows}
              unit="MB"
              axisLabel="Max peak RSS, lower is better"
              scale="linear"
            />
          </article>

          <article className="benchmark-card">
            <div className="chart-copy">
              <span className="card-label">Count agreement</span>
              <h3>Workflow comparisons are not treated as oracles.</h3>
              <p>
                MAGeCK and guide-counter rows are useful adoption comparisons.
                The assignment oracle remains the native exhaustive scan.
              </p>
            </div>
            <AgreementChart rows={agreementRows} />
          </article>

          <article className="benchmark-card">
            <div className="chart-copy">
              <span className="card-label">Multi-sample CRISPR</span>
              <h3>crispr-count keeps one assignment per read.</h3>
              <p>
                The current scaling artifact covers 2, 4, and 8 Yusa samples
                with threaded DotMatch Hamming k=1. DotMatch records zero
                overcount reads; guide-counter's multi-offset lane reports more
                counted assignments than input reads.
              </p>
            </div>
            <HorizontalBarChart
              rows={scalingRows}
              unit="reads"
              axisLabel="Overcount reads, lower is better"
              scale="linear"
            />
          </article>
        </div>
      </section>

      <section id="use-cases" className="section use-cases">
        <div className="section-heading">
          <h2>One engine, several short-DNA workflows.</h2>
          <p>
            These are the current implemented surfaces from the CLI and
            benchmark docs. Genome mapping, SAM/BAM, CIGAR, wildcard N
            semantics, CBCL/NovaSeq input, and production wheels are not
            current claims.
          </p>
        </div>
        <div className="usecase-grid">
          {useCases.map((item) => (
            <article key={item.title} className="usecase">
              <span className="usecase-dot" />
              <h3>{item.title}</h3>
              <p>{item.body}</p>
            </article>
          ))}
        </div>
      </section>

      <section id="comparison" className="section comparison">
        <div className="section-heading">
          <h2>Built for known-target reads.</h2>
          <p>
            Give DotMatch a target library and sequencing reads. It assigns
            reads, reports ambiguity instead of hiding it, and emits the files
            needed to reproduce, audit, validate, and debug the result.
          </p>
        </div>
        <div className="comparison-layout">
          <div className="comparison-table" role="table" aria-label="DotMatch current product surface">
            <div role="row" className="table-head">
              <span>Surface</span>
              <span>Current support</span>
              <span>Boundary</span>
            </div>
            <Row
              a="Core assignment"
              b="Exact global edit distance, threshold queries, k=0/k=1 index, best or radius ambiguity policy"
              c="No semi-global/infix alignment, traceback, CIGAR, or wildcard N semantics yet"
            />
            <Row
              a="CRISPR counting"
              b="FASTQ/FASTQ.gz count, crispr-count wrapper, MAGeCK matrix, QC JSON/TSV, HTML report"
              c="Claim is Yusa-style guide counting, not broad screen analysis"
            />
            <Row
              a="Inline demux"
              b="Fixed-position single-end FASTQ/FASTQ.gz barcode splitting with ambiguous and unmatched outputs"
              c="Synthetic fixture is a smoke benchmark; real public barcode rows are still gated"
            />
            <Row
              a="Classic BCL demux"
              b="RunInfo, SampleSheet v1/v2 data sections, classic per-cycle BCL(.gz), filter files, sample FASTQ.gz"
              c="CBCL/NovaSeq-style input remains a planned milestone"
            />
            <Row
              a="Trust layer"
              b="audit, validate, optional Edlib oracle, inspect-unmatched, public TSV/JSON schemas"
              c="Production workbench and hosted infrastructure are future proprietary layers"
            />
          </div>
          <div className="evidence-notes">
            <article>
              <span>Synthetic inline demux</span>
              <strong>918k reads/s</strong>
              <p>
                Current smoke fixture: 20k reads, 4 barcodes, Hamming k=1,
                matching Cutadapt assigned/unmatched totals. This is pipeline
                evidence, not a state-of-the-art barcode claim.
              </p>
            </article>
            <article>
              <span>Classic BCL demux</span>
              <strong>61k clusters/s</strong>
              <p>
                Public 10x tiny-BCL row: 2.14M clusters, 132 cycles, one sample.
                bcl2fastq comparator row records 20.8k clusters/s and zero
                validation mismatches; broader BCL claims remain gated.
              </p>
            </article>
          </div>
        </div>
      </section>

      <section id="workflow" className="section workflow">
        <div className="workflow-copy">
          <h2>Designed for reproducible command-line workflows.</h2>
          <p>
            The C core ships as a CLI, static/shared library, and Python ctypes
            bindings. The public workflow emits count matrices, FASTQ splits,
            QC files, assignment diagnostics, audit tables, validation summaries,
            and deterministic HTML reports.
          </p>
        </div>
        <div className="terminal" aria-label="DotMatch commands">
          <div className="terminal-bar">
            <span />
            <span />
            <span />
          </div>
          {commands.map((command) => (
            <code key={command}>
              <span>$</span> {command}
            </code>
          ))}
        </div>
      </section>

      <section className="section final-cta">
        <h2>Known-target assignment without hiding the hard cases.</h2>
        <p>
          Use DotMatch when short reads must be assigned to a known target set,
          one-edit correction should be exact, and ambiguous, unsafe, or
          unmatched cases need to stay visible.
        </p>
        <a className="button primary" href="#benchmarks">
          Review the evidence
        </a>
      </section>
    </main>
  );
}

function Row({ a, b, c }: { a: string; b: string; c: string }) {
  return (
    <div role="row">
      <span data-label="Surface">{a}</span>
      <span data-label="Current support">{b}</span>
      <span data-label="Boundary">{c}</span>
    </div>
  );
}

function HorizontalBarChart({
  rows,
  unit,
  axisLabel,
  scale
}: {
  rows: readonly { label: string; value: number; tone: string }[];
  unit: string;
  axisLabel: string;
  scale: "linear" | "log";
}) {
  const max = Math.max(...rows.map((row) => row.value));
  const logFloor = 1;
  const ticks = scale === "log" ? logTicks(max) : [0, max * 0.25, max * 0.5, max * 0.75, max];
  const ariaSummary = rows
    .map((row) => `${row.label}: ${formatNumber(row.value)} ${unit}`)
    .join("; ");

  function width(value: number) {
    if (scale === "log") {
      const min = Math.log10(logFloor);
      const range = Math.log10(max) - min || 1;
      return ((Math.log10(Math.max(value, logFloor)) - min) / range) * 100;
    }

    return (value / max) * 100;
  }

  return (
    <div className="native-chart" role="img" aria-label={`${axisLabel}. ${ariaSummary}.`}>
      <div className="chart-axis-label">{axisLabel}</div>
      <div className="chart-plot">
        <div className="chart-gridlines" aria-hidden="true">
          {ticks.map((tick) => {
            const left = scale === "log" ? width(tick) : (tick / max) * 100;
            return <span key={tick} style={{ left: `${Math.min(left, 100)}%` }} />;
          })}
        </div>
        <div className="bar-list">
          {rows.map((row) => (
            <div className="bar-row" key={row.label}>
              <div className="bar-meta">
                <span>{row.label}</span>
                <strong>
                  {formatNumber(row.value)}
                  <em>{unit}</em>
                </strong>
              </div>
              <div className="bar-track">
                <span
                  className={`tone-${row.tone}`}
                  style={{ width: `${Math.max(width(row.value), 1.5)}%` }}
                />
              </div>
            </div>
          ))}
        </div>
        <div className="chart-ticks" aria-hidden="true">
          {ticks.map((tick) => {
            const left = scale === "log" ? width(tick) : (tick / max) * 100;
            return (
              <span key={tick} style={{ left: `${Math.min(left, 100)}%` }}>
                {formatCompact(tick)}
              </span>
            );
          })}
        </div>
      </div>
    </div>
  );
}

function AgreementChart({
  rows
}: {
  rows: readonly { label: string; value: number; tone: string }[];
}) {
  const ticks = [0, 0.25, 0.5, 0.75, 1];
  const ariaSummary = rows
    .map((row) => `${row.label}: Pearson ${row.value.toFixed(3)}`)
    .join("; ");

  return (
    <div
      className="native-chart agreement-chart"
      role="img"
      aria-label={`Pearson agreement by workflow. ${ariaSummary}.`}
    >
      <div className="chart-axis-label">Pearson correlation by guide count table</div>
      <div className="chart-plot">
        <div className="chart-gridlines" aria-hidden="true">
          {ticks.map((tick) => (
            <span key={tick} style={{ left: `${tick * 100}%` }} />
          ))}
        </div>
        <div className="agreement-list">
          {rows.map((row) => (
            <div className="agreement-row" key={row.label}>
              <div className="agreement-meta">
                <span>{row.label}</span>
                <strong>{row.value.toFixed(3)}</strong>
              </div>
              <div className="bar-track">
                <span className={`tone-${row.tone}`} style={{ width: `${row.value * 100}%` }} />
              </div>
            </div>
          ))}
        </div>
        <div className="chart-ticks" aria-hidden="true">
          {ticks.map((tick) => (
            <span key={tick} style={{ left: `${tick * 100}%` }}>
              {tick.toFixed(tick === 0 || tick === 1 ? 0 : 2)}
            </span>
          ))}
        </div>
      </div>
    </div>
  );
}

function formatNumber(value: number) {
  return value.toLocaleString(undefined, {
    maximumFractionDigits: value < 100 ? 1 : 0
  });
}

function formatCompact(value: number) {
  if (value >= 1000000) {
    const scaled = value / 1000000;
    return `${scaled >= 10 ? Math.round(scaled) : Number(scaled.toFixed(1))}M`;
  }

  if (value >= 1000) {
    const scaled = value / 1000;
    return `${scaled >= 10 ? Math.round(scaled) : Number(scaled.toFixed(1))}k`;
  }

  return value.toLocaleString(undefined, {
    maximumFractionDigits: value < 10 ? 1 : 0
  });
}

function logTicks(max: number) {
  const ticks = [];
  const topPower = Math.ceil(Math.log10(Math.max(max, 1)));

  for (let power = 0; power <= topPower; power += 1) {
    ticks.push(10 ** power);
  }

  return ticks;
}

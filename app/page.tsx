const repoUrl = "https://github.com/dnncha/dotmatch";
const citationUrl = `${repoUrl}/blob/main/CITATION.cff`;
const methodsUrl = `${repoUrl}/blob/main/docs/methods-and-citation.md`;
const publicCrisprUrl = `${repoUrl}/blob/main/docs/benchmarks/public_crispr/README.md`;

const proof = [
  ["Best fit", "fixed target lists", "Guides, barcodes, primers, adapters, panels, and whitelist-style sequences where the candidates are already known."],
  ["Correctness rule", "index matches scan", "The fast path is tested against the native exhaustive scan for the same targets, error allowance, and ambiguity policy."],
  ["Public CRISPR data", "5 repeated rows", "The checked-in Yusa/MAGeCK rows cover 10k and 100k reads per sample, with count agreement and independent edit-distance checks."],
  ["Repository contents", "C, CLI, Python", "Core code, bindings, tests, scripts, reports, schemas, and raw benchmark tables live in the repo."]
];

const decisionCards = [
  {
    title: "Use DotMatch when you have",
    items: [
      "CRISPR guide-counting FASTQs",
      "inline barcode reads",
      "known primer, panel, or whitelist targets",
      "classic per-cycle BCL demultiplexing jobs"
    ]
  },
  {
    title: "DotMatch gives you",
    items: [
      "one assignment per read",
      "explicit ambiguous and unmatched reads",
      "one-base mismatch or indel rescue",
      "MAGeCK-compatible count matrices and QC tables"
    ]
  },
  {
    title: "Do not use DotMatch for",
    items: [
      "genome alignment or variant calling",
      "SAM/BAM/CIGAR output",
      "downstream CRISPR screen statistics",
      "CBCL/NovaSeq demultiplexing or wildcard N semantics"
    ]
  }
];

const translations = [
  ["known targets", "a fixed guide, barcode, primer, whitelist, or panel sequence list"],
  ["Hamming k=1", "allow one mismatch, no indels"],
  ["Levenshtein k=1", "allow one substitution, insertion, or deletion"],
  ["ambiguous", "reads that match multiple targets are reported, not forced into a guide or barcode"],
  ["peak RSS", "peak memory use"],
  ["Edlib validation", "checked against an independent edit-distance implementation"]
];

const audienceCards = [
  {
    title: "CRISPR screen users",
    body: "Count guides from FASTQ/FASTQ.gz into MAGeCK-compatible matrices, with exact, rescued, ambiguous, and unmatched reads in the QC."
  },
  {
    title: "Sequencing cores",
    body: "Demultiplex fixed-position inline barcodes while keeping ambiguous and unmatched reads available for review."
  },
  {
    title: "Bioinformatics developers",
    body: "Use the C core, CLI, Python bindings, schemas, validation commands, and raw benchmark artifacts."
  },
  {
    title: "Methods reviewers",
    body: "Inspect the claim gates, raw CSVs, exact commands, and validation against exhaustive or Edlib checks."
  }
];

const workflowStatusRows = [
  ["CRISPR guide counting", "Validated now", "Public MAGeCK/Yusa repeated rows, count agreement, Edlib validation, and raw command tables."],
  ["Target-library audit", "Supported", "CLI tests, schemas, and validation commands for unsafe one-edit libraries."],
  ["Inline barcode demux", "Supported, bounded", "Fixed-position FASTQ demux and smoke benchmark; broader real-data comparison still gated."],
  ["Classic BCL demux", "Milestone", "Public 10x tiny-BCL row; CBCL/NovaSeq input is not supported yet."],
  ["Genome alignment", "Not supported", "Out of scope: no SAM/BAM/CIGAR, reference mapping, or variant calling."]
];

const workflowChoiceRows = [
  ["Downstream CRISPR screen statistics", "MAGeCK or another downstream screen-analysis tool"],
  ["FASTQ-to-guide count matrix with explicit ambiguity QC", "DotMatch"],
  ["Genome or transcriptome reference mapping", "Bowtie2, BWA, minimap-style tools, not DotMatch"],
  ["Adapter trimming", "Cutadapt-style tools, not DotMatch"],
  ["Known short target assignment with exact one-edit semantics", "DotMatch"]
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

const basePath = process.env.NEXT_PUBLIC_BASE_PATH ?? "";
const heroWorkflowImage = `${basePath}/dotmatch-hero-workflow.png`;

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
          <a href="#install">Install</a>
          <a href="#use-cases">Use cases</a>
          <a href="#cite">Cite</a>
          <a href={repoUrl}>GitHub</a>
        </nav>
        <a className="header-cta" href={repoUrl}>Source</a>
      </header>

      <section id="top" className="hero">
        <div className="hero-copy">
          <h1>DotMatch</h1>
          <p className="hero-lede">
            Auditable assignment of short sequencing reads to known DNA targets.
          </p>
          <p className="hero-text">
            DotMatch counts CRISPR guides, assigns inline barcodes, and matches
            short primer, panel, or whitelist reads against a known target list.
            It supports exact matching, one-mismatch rescue, and one-base
            insertion or deletion rescue while reporting ambiguous reads instead
            of silently guessing.
          </p>
          <p className="hero-note">
            <strong>Best-supported today: CRISPR guide counting</strong> from
            public MAGeCK/Yusa FASTQs, with MAGeCK-compatible count matrices, QC
            summaries, raw benchmark artifacts, and validation checks.
          </p>
          <div className="hero-actions">
            <a href="#benchmarks" className="button primary">
              Read the CRISPR evidence
            </a>
            <a href="#install" className="button secondary">
              Install
            </a>
            <a href={repoUrl} className="button secondary">
              GitHub
            </a>
          </div>
        </div>
        <div className="hero-panel" aria-label="DotMatch benchmark summary">
          <div className="panel-topline">
            <span>v0.1.0</span>
            <span>CRISPR-led evidence</span>
          </div>
          <figure className="hero-visual">
            <img
              src={heroWorkflowImage}
              alt="Short DNA reads flowing into known target assignments while ambiguous and unmatched reads stay visible"
              decoding="async"
              fetchPriority="high"
            />
            <figcaption>
              Reads move into known targets; ambiguous and unmatched lanes stay visible.
            </figcaption>
          </figure>
          <div className="metric-grid">
            <div>
              <strong>87,437</strong>
              <span>public MAGeCK/Yusa guides in the current benchmark fixture</span>
            </div>
            <div>
              <strong>0</strong>
              <span>Edlib validation mismatches across 2,000 checked reads</span>
            </div>
            <div>
              <strong>331k</strong>
              <span>reads/s for one-mismatch CRISPR guide counting</span>
            </div>
            <div>
              <strong>28.7 MB</strong>
              <span>peak memory use in the repeated public CRISPR rows</span>
            </div>
          </div>
          <div className="sequence-rail" aria-hidden="true">
            {Array.from({ length: 64 }).map((_, i) => (
              <span key={i} className={i % 7 === 0 ? "hot" : i % 5 === 0 ? "cool" : ""} />
            ))}
          </div>
        </div>
      </section>

      <section className="section decision-section" aria-label="DotMatch adoption guide">
        <div className="section-heading">
          <h2>Use it when hidden assignment choices matter.</h2>
          <p>
            Most DotMatch jobs start as FASTQ reads and a target table. The
            point is not only speed; it is making corrected, ambiguous, and
            unmatched reads visible enough to audit.
          </p>
        </div>
        <div className="decision-grid">
          {decisionCards.map((card) => (
            <article key={card.title} className="decision-card">
              <h3>{card.title}</h3>
              <ul>
                {card.items.map((item) => (
                  <li key={item}>{item}</li>
                ))}
              </ul>
            </article>
          ))}
        </div>
        <div className="translation-grid" aria-label="Biology translations for DotMatch terms">
          {translations.map(([term, meaning]) => (
            <div key={term}>
              <span>{term}</span>
              <p>{meaning}</p>
            </div>
          ))}
        </div>
      </section>

      <section className="section example-section">
        <div className="section-heading">
          <h2>One CRISPR run, from FASTQ to QC.</h2>
          <p>
            This is the practical shape of the best-supported workflow: reads
            in, a guide-by-sample count matrix out, and a small set of QC files
            that say what happened to every assignment class.
          </p>
        </div>
        <div className="example-layout">
          <article className="example-card">
            <span className="card-label">Minimal example</span>
            <pre><code>{`dotmatch crispr-count \\
  --library yusa_library.csv \\
  --samples samples.tsv \\
  --guide-start 23 \\
  --guide-length 19 \\
  --k 1 \\
  --metric levenshtein \\
  --indel-window 1 \\
  --out counts.mageck.tsv \\
  --summary qc.json \\
  --report report.html`}</code></pre>
            <div className="output-list" aria-label="DotMatch CRISPR outputs">
              <code>counts.mageck.tsv</code>
              <span>guide x sample count matrix</span>
              <code>qc.json</code>
              <span>exact, rescued, ambiguous, and unmatched reads</span>
              <code>report.html</code>
              <span>archived run report</span>
            </div>
          </article>
          <article className="ambiguity-example">
            <span className="card-label">Why ambiguity is explicit</span>
            <pre><code>{`Read:    ACGTACGT
Guide A: ACGTACGA   distance 1
Guide B: ACGTACGC   distance 1

Some tools may pick or double-count.
DotMatch reports: ambiguous`}</code></pre>
            <p>
              Ambiguous reads are not silently counted into a guide or barcode.
              They stay available for QC and diagnosis.
            </p>
          </article>
        </div>
      </section>

      <section id="install" className="section launch-section">
        <div className="section-heading">
          <h2>Start from the repo. Cite the exact release.</h2>
          <p>
            DotMatch is not on every package channel yet, so the honest first
            install is from source. Current distribution: source install and
            release artifacts. Coming next: PyPI, Bioconda, Docker/Singularity,
            Zenodo DOI.
          </p>
        </div>
        <div className="launch-grid">
          <article className="launch-card">
            <span className="card-label">Build it locally</span>
            <h3>Clone the repo and run the release check.</h3>
            <pre><code>{`git clone https://github.com/dnncha/dotmatch.git
cd dotmatch
make
python3 -m pip install .
make repository-ready`}</code></pre>
            <a href={repoUrl}>Open GitHub</a>
          </article>

          <article id="cite" className="launch-card">
            <span className="card-label">Cite it</span>
            <h3>Use the release citation and a matching methods sentence.</h3>
            <p>
              If DotMatch helps an analysis, cite the software release. The
              methods note has short wording for CRISPR guide counting,
              one-edit Levenshtein rescue, and Hamming-only comparisons.
            </p>
            <div className="link-stack">
              <a href={citationUrl}>CITATION.cff</a>
              <a href={methodsUrl}>Methods and citation notes</a>
            </div>
          </article>

          <article className="launch-card">
            <span className="card-label">Check the data</span>
            <h3>The main public comparison is deliberately narrow.</h3>
            <p>
              The public CRISPR benchmark is the best-supported comparison
              today: Yusa-style guide counting, checked-in rows, and validation
              against the assignment oracle.
            </p>
            <div className="link-stack">
              <a href={publicCrisprUrl}>Public CRISPR benchmark report</a>
              <a href="#benchmarks">Review benchmark summary</a>
            </div>
          </article>
        </div>
      </section>

      <section id="benchmarks" className="section proof-section">
        <div className="section-heading">
          <h2>The public CRISPR evidence in plain English.</h2>
          <p>
            We are keeping the claims narrow for v0.1.0. On repeated public
            MAGeCK/Yusa CRISPR guide-counting rows, DotMatch Hamming k=1
            processed about 331k reads/s using about 28.7 MB peak memory;
            guide-counter processed about 195k reads/s using about 529 MB, and
            MAGeCK exact count processed about 93k reads/s using about 159 MB.
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
              <span className="card-label">Public CRISPR benchmark</span>
              <h3>The Yusa rows are in the repo.</h3>
              <p>
                These rows are not a leaderboard. They are the first public case
                we can rerun and inspect: five 100k-record/sample repeats for
                DotMatch, MAGeCK, and guide-counter, with exact, Hamming, and
                Levenshtein kept separate. Edlib validation checks 2,000 reads
                with zero mismatches against an independent edit-distance
                implementation.
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
              <h3>k=1 Levenshtein usually checks only a few candidates.</h3>
              <p>
                On the public Yusa rows, the index sends about 2.822 candidate
                targets per read to exact verification, out of an 87,437-guide
                library. In biology terms, that lane allows one substitution,
                insertion, or deletion.
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
              <h3>The CRISPR counter stays small.</h3>
              <p>
                The repeated Yusa runs put DotMatch Hamming and exact lanes
                around 28.7 MB peak memory use. guide-counter is around 528.7
                MB on the same fixture.
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
              <h3>Comparator counts are useful, but not oracles.</h3>
              <p>
                MAGeCK and guide-counter help us compare familiar workflows.
                Correctness is checked against exhaustive assignment and Edlib,
                not whichever external tool happens to agree.
              </p>
            </div>
            <AgreementChart rows={agreementRows} />
          </article>

          <article className="benchmark-card">
            <div className="chart-copy">
              <span className="card-label">Multi-sample CRISPR</span>
              <h3>One read, one counted assignment.</h3>
              <p>
                In the 2/4/8-sample scaling run, DotMatch records zero
                overcount reads. The guide-counter multi-offset lane can count
                more assignments than reads.
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
          <h2>Who this is for.</h2>
          <p>
            The same engine serves a few different readers. The strongest
            adoption path today is CRISPR guide counting, but the audit trail is
            useful anywhere short reads must land on a fixed target list.
          </p>
        </div>
        <div className="usecase-grid">
          {audienceCards.map((item) => (
            <article key={item.title} className="usecase">
              <span className="usecase-dot" />
              <h3>{item.title}</h3>
              <p>{item.body}</p>
            </article>
          ))}
        </div>
      </section>

      <section className="section status-section">
        <div className="section-heading">
          <h2>Validated now, available but early, or out of scope.</h2>
          <p>
            The project has several surfaces. This table keeps the strongest
            public claim separate from smoke tests and future distribution work.
          </p>
        </div>
        <div className="status-table" role="table" aria-label="DotMatch workflow maturity">
          <div role="row" className="table-head">
            <span>Workflow</span>
            <span>Status</span>
            <span>Evidence level</span>
          </div>
          {workflowStatusRows.map(([workflow, status, evidence]) => (
            <div role="row" key={workflow}>
              <span data-label="Workflow">{workflow}</span>
              <span data-label="Status">{status}</span>
              <span data-label="Evidence level">{evidence}</span>
            </div>
          ))}
        </div>
      </section>

      <section id="comparison" className="section comparison">
        <div className="section-heading">
          <h2>How to place it beside familiar tools.</h2>
          <p>
            DotMatch is not trying to replace every CRISPR, trimming, or mapping
            tool. It is the assignment and QC step when a read should map to one
            known short DNA target.
          </p>
        </div>
        <div className="comparison-layout">
          <div className="comparison-table" role="table" aria-label="DotMatch current CLI support">
            <div role="row" className="table-head">
              <span>Need</span>
              <span>Use</span>
            </div>
            {workflowChoiceRows.map(([need, tool]) => (
              <div role="row" key={need}>
                <span data-label="Need">{need}</span>
                <span data-label="Use">{tool}</span>
              </div>
            ))}
          </div>
          <div className="evidence-notes">
            <article>
              <span>One read policy</span>
              <strong>0</strong>
              <p>
                DotMatch records zero overcount reads in the 2/4/8-sample CRISPR
                scaling run. Ambiguity is reported instead of being hidden in a
                count total.
              </p>
            </article>
            <article>
              <span>Bounded demux evidence</span>
              <strong>918k</strong>
              <p>
                The inline demux row is a 20k-read smoke fixture with four
                barcodes and Hamming k=1. It is useful pipeline evidence, not a
                broad barcode benchmark.
              </p>
            </article>
          </div>
        </div>
      </section>

      <section id="workflow" className="section workflow">
        <div className="workflow-copy">
          <h2>Command-line first.</h2>
          <p>
            DotMatch is a small C/Python tool with a CLI and Python ctypes
            bindings. Runs can write count matrices, FASTQ splits, QC tables,
            assignment diagnostics, audit files, validation summaries, and
            self-contained HTML reports.
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
        <h2>For short reads with known targets and real QC stakes.</h2>
        <p>
          Use DotMatch when exact one-edit assignment matters, when ambiguous or
          unmatched reads are as important as the counts, and when another lab
          should be able to inspect how the calls were made.
        </p>
        <a className="button primary" href="#benchmarks">
          Review the evidence
        </a>
      </section>
    </main>
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

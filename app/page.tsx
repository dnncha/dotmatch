const repoUrl = "https://github.com/dnncha/dotmatch";
const citationUrl = `${repoUrl}/blob/main/CITATION.cff`;
const methodsUrl = `${repoUrl}/blob/main/docs/methods-and-citation.md`;
const packagingUrl = `${repoUrl}/blob/main/docs/packaging.md`;
const publicCrisprUrl = `${repoUrl}/blob/main/docs/benchmarks/public_crispr/README.md`;
const barcodeScienceUrl = `${repoUrl}/blob/main/docs/barcode-science-readiness.md`;
const barcodeBenchmarkUrl = `${repoUrl}/blob/main/docs/benchmarks/barcode_demux/README.md`;
const biocondaPrUrl = "https://github.com/bioconda/bioconda-recipes/pull/65367";

const proof = [
  ["4 outcomes", "unique, ambiguous, none, invalid", "Every read lands in an auditable assignment class."],
  ["8 diagnoses", "barcode autopsy", "Wrong offset, unsafe rescue, collisions, unmatched classes, and more."],
  ["5 lanes", "public evidence", "Fixed-window barcode/CRISPR checks are documented with scoped claims."],
  ["0 silent picks", "ambiguity retained", "Ambiguous reads stay visible instead of being guessed into a sample."]
];

const decisionCards = [
  {
    title: "Use DotMatch when you have",
    items: [
      "fixed-window barcode FASTQs",
      "CRISPR guide-counting reads",
      "known primer, panel, or whitelist targets",
      "feature-barcode or amplicon slices"
    ]
  },
  {
    title: "DotMatch gives you",
    items: [
      "one assignment per read",
      "explicit ambiguous and unmatched reads",
      "unsafe one-edit correction warnings",
      "HTML, TSV, JSON, and workflow artifacts"
    ]
  },
  {
    title: "Do not use DotMatch for",
    items: [
      "genome alignment or variant calling",
      "SAM/BAM/CIGAR output",
      "downstream CRISPR screen statistics",
      "BCL Convert replacement workflows"
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
    title: "Barcode assay owners",
    body: "Find the barcode window, audit the barcode library, demultiplex fixed-position reads, and explain unassigned reads."
  },
  {
    title: "Sequencing cores",
    body: "Turn undetermined inline-barcode lanes into reports that show wrong offsets, collisions, and unsafe rescue choices."
  },
  {
    title: "CRISPR screen users",
    body: "Count guides from FASTQ/FASTQ.gz into MAGeCK-compatible matrices, with exact, rescued, ambiguous, and unmatched reads in the QC."
  },
  {
    title: "Methods reviewers",
    body: "Inspect claim gates, raw CSVs, exact commands, deterministic reports, and validation against exhaustive or Edlib checks."
  }
];

const workflowStatusRows = [
  ["Barcode autopsy", "Supported in the repository", "Offset scans, safety audit, unmatched summaries, provenance, and conservative public-surface checks."],
  ["Inline barcode demux", "Comparator-backed, bounded", "Fixed-position exact-prefix public lane with documented Cutadapt/hash-splitter semantics."],
  ["CRISPR guide counting", "Validated now", "Public MAGeCK/Yusa repeated rows, count agreement, Edlib validation, and raw command tables."],
  ["Target-library audit", "Supported", "CLI tests, schemas, and validation commands for unsafe one-edit libraries."],
  ["Classic BCL demux", "Early milestone", "Public 10x tiny-BCL row; BCL Convert and CBCL/NovaSeq inputs remain outside the flagship story."],
  ["Genome alignment", "Not supported", "No SAM/BAM/CIGAR, reference mapping, or variant calling."]
];

const workflowChoiceRows = [
  ["Downstream CRISPR screen statistics", "MAGeCK or another downstream screen-analysis tool"],
  ["FASTQ-to-guide count matrix with explicit ambiguity QC", "DotMatch"],
  ["Genome or transcriptome reference mapping", "Bowtie2, BWA, minimap-style tools, not DotMatch"],
  ["Adapter trimming", "Cutadapt-style tools, not DotMatch"],
  ["Known short target assignment with exact one-edit semantics", "DotMatch"]
];

const evidenceNotes = [
  ["Correctness rule", "index matches scan", "The fast path is tested against the native exhaustive scan for the same targets, error allowance, and ambiguity policy."],
  ["Best fit", "fixed target lists", "Guides, barcodes, primers, adapters, panels, and whitelist-style sequences where the candidates are already known."],
  ["Repository contents", "C, CLI, Python", "Core code, bindings, tests, scripts, reports, schemas, and raw benchmark tables live in the repo."]
];

const commands = [
  "dotmatch barcode autopsy --barcodes barcodes.tsv --reads pooled.fastq.gz --scan-starts 0:12 --k-values 0,1 --out-dir autopsy",
  "dotmatch barcode infer --barcodes barcodes.tsv --reads pooled.fastq.gz --scan-starts 0:30 --sample-reads 100000 --out inference.tsv",
  "dotmatch barcode demux --barcodes barcodes.tsv --reads pooled.fastq.gz --barcode-start 1 --barcode-length auto --k 1 --metric hamming --max-correction-qual 20 --out-dir demuxed --report report.html",
  "dotmatch crispr-count --library guides.csv --samples samples.tsv --guide-start 23 --guide-length 19 --k 1 --metric levenshtein --indel-window 1 --out counts.mageck.tsv --summary qc.json",
  "dotmatch assay run assay.toml",
  "dotmatch validate --targets guides.tsv --reads sample.fastq.gz --target-start 23 --target-length 19 --k 1 --indel-window 1 --oracle edlib --sample 100000"
];

const autopsyArtifacts = [
  ["report.html", "scientist-readable barcode autopsy report"],
  ["findings.tsv", "wrong offset, unsafe rescue, and collision findings"],
  ["offset_scan.tsv", "candidate windows ranked by assignment evidence"],
  ["correction_safety.tsv", "whether one-edit rescue is safe for this barcode set"],
  ["top_unmatched.tsv", "high-count unassigned barcode sequences"],
  ["provenance.json", "commands, versions, thresholds, and artifacts"]
];

const autopsyFindings = [
  ["wrong offset", "Detects a likely leading base, primer scar, or shifted barcode window."],
  ["unsafe correction", "Shows barcode pairs or clusters that make k=1 rescue scientifically risky."],
  ["ambiguous collision", "Keeps reads that match multiple barcodes out of forced assignments."],
  ["unmatched classes", "Separates low-complexity, distant, reverse-complement, and quality-gated failures."]
];

const reportPreviewRows = [
  ["unique", "assigned to exactly one barcode or guide", "counted or split"],
  ["ambiguous", "compatible with multiple targets", "reported, not forced"],
  ["none", "outside the configured edit radius", "sent to unmatched diagnostics"],
  ["invalid", "window could not be extracted", "kept visible in QC"]
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

const basePath = process.env.NEXT_PUBLIC_BASE_PATH ?? "";
const barcodeAutopsyImage = `${basePath}/dotmatch-barcode-autopsy.png`;

export default function Home() {
  return (
    <main>
      <header className="site-header">
        <a className="brand" href="#top" aria-label="DotMatch home">
          <span className="brand-mark" />
          DotMatch
        </a>
        <nav aria-label="Primary navigation">
          <a href="#autopsy">Autopsy</a>
          <a href="#benchmarks">Benchmarks</a>
          <a href="#use-cases">Use cases</a>
          <a href="#install">Install</a>
          <a href="#cite">Cite</a>
          <a href={repoUrl}>GitHub</a>
        </nav>
        <a className="header-cta" href={repoUrl}>Source</a>
      </header>

      <section id="top" className="hero">
        <div className="hero-copy">
          <h1>DotMatch</h1>
          <p className="hero-lede">
            Barcode autopsy for fixed-window FASTQ assays.
          </p>
          <p className="hero-text">
            DotMatch turns barcode, guide, feature-tag, primer, and whitelist
            assignment into an auditable report: what matched, what was
            ambiguous, what failed by offset, what is unsafe to rescue, and
            exactly which commands produced the result.
          </p>
          <p className="hero-note">
            <strong>Flagship story: post-FASTQ fixed-window assignment.</strong>{" "}
            DotMatch does not replace BCL Convert, genome aligners, or general
            adapter trimming. It makes known-target assignment explainable.
          </p>
          <div className="hero-actions">
            <a href="#autopsy" className="button primary">
              See Barcode Autopsy
            </a>
            <a href={barcodeScienceUrl} className="button secondary">
              Run the science gate
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
            <span>v0.1.1</span>
            <span>fixed-window evidence</span>
          </div>
          <figure className="hero-visual">
            <img
              src={barcodeAutopsyImage}
              alt="DotMatch barcode autopsy workflow from FASTQ through fixed window audit into unique, ambiguous, none, invalid, report, and provenance outputs"
              decoding="async"
              fetchPriority="high"
            />
            <figcaption>
              FASTQ reads become unique, ambiguous, none, and invalid outcomes,
              with reports and provenance kept beside the split outputs.
            </figcaption>
          </figure>
          <div className="metric-grid">
            <div>
              <strong>4</strong>
              <span>assignment outcomes, including ambiguous and invalid reads</span>
            </div>
            <div>
              <strong>8</strong>
              <span>barcode failure classes covered by the fixture catalog</span>
            </div>
            <div>
              <strong>1.37M</strong>
              <span>reads/s on the bounded exact-prefix barcode lane</span>
            </div>
            <div>
              <strong>0</strong>
              <span>silent ambiguity when one-edit correction is unsafe</span>
            </div>
          </div>
          <div className="sequence-rail" aria-hidden="true">
            {Array.from({ length: 64 }).map((_, i) => (
              <span key={i} className={i % 7 === 0 ? "hot" : i % 5 === 0 ? "cool" : ""} />
            ))}
          </div>
        </div>
      </section>

      <section className="evidence-strip" aria-label="DotMatch public CRISPR evidence summary">
        {proof.map(([label, value, detail]) => (
          <article key={label}>
            <strong>{label}</strong>
            <span>{value}</span>
            <p>{detail}</p>
          </article>
        ))}
      </section>

      <section id="autopsy" className="section autopsy-section">
        <div className="section-heading">
          <h2>The demo users understand in one run.</h2>
          <p>
            Barcode Autopsy answers the question demultiplexing logs usually
            leave open: why were these reads not assigned, and is correction
            safe for this barcode set? Speed is shown only after the comparator
            semantics are documented.
          </p>
        </div>
        <div className="autopsy-layout">
          <article className="autopsy-command">
            <span className="card-label">Flagship command</span>
            <pre><code>{`dotmatch barcode autopsy \\
  --barcodes barcodes.tsv \\
  --reads pooled.fastq.gz \\
  --scan-starts 0:12 \\
  --k-values 0,1 \\
  --out-dir autopsy`}</code></pre>
            <p>
              The command scans plausible windows, audits the barcode library,
              writes split-read evidence, and produces report files that can be
              attached to a sequencing handoff or workflow run.
            </p>
          </article>
          <div className="artifact-grid" aria-label="Barcode autopsy outputs">
            {autopsyArtifacts.map(([name, detail]) => (
              <article key={name}>
                <code>{name}</code>
                <p>{detail}</p>
              </article>
            ))}
          </div>
        </div>
        <div className="finding-list" aria-label="Barcode autopsy diagnosis examples">
          {autopsyFindings.map(([label, detail]) => (
            <article key={label}>
              <span>{label}</span>
              <p>{detail}</p>
            </article>
          ))}
        </div>
      </section>

      <section className="section report-section">
        <div className="section-heading">
          <h2>Readable for the bench. Exact for the workflow.</h2>
          <p>
            The primary artifact is an HTML report that explains the run in
            plain language, backed by TSV and JSON files that a pipeline can
            archive, compare, and summarize.
          </p>
        </div>
        <div className="report-preview" aria-label="DotMatch report outcome preview">
          <div className="report-copy">
            <h3>Every read keeps its assignment reason.</h3>
            <p>
              DotMatch does not turn a failed demultiplexing run into one
              undifferentiated “undetermined” bucket. It separates ambiguous
              rescue, wrong windows, invalid slices, and true no-match reads so
              a scientist can decide whether to change the assay spec or rerun
              the sample.
            </p>
          </div>
          <div className="report-table" role="table" aria-label="Assignment outcome meanings">
            <div role="row" className="table-head">
              <span>Outcome</span>
              <span>Meaning</span>
              <span>Action</span>
            </div>
            {reportPreviewRows.map(([outcome, meaning, action]) => (
              <div role="row" key={outcome}>
                <span data-label="Outcome"><code>{outcome}</code></span>
                <span data-label="Meaning">{meaning}</span>
                <span data-label="Action">{action}</span>
              </div>
            ))}
          </div>
        </div>
      </section>

      <section id="benchmarks" className="section proof-section">
        <div className="section-heading">
          <h2>The public evidence in plain English.</h2>
          <p>
            We are keeping the claims narrow for v0.1.1. The barcode lane
            documents fixed-position exact-prefix semantics against Cutadapt and
            a hash splitter, while the CRISPR rows document guide-counting
            behavior. On repeated public
            MAGeCK/Yusa CRISPR guide-counting rows, DotMatch Hamming k=1
            processed about 331k reads/s using about 28.7 MB peak memory;
            guide-counter processed about 195k reads/s using about 529 MB, and
            MAGeCK exact count processed about 93k reads/s using about 159 MB.
          </p>
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
            <div className="link-stack compact">
              <a href={barcodeBenchmarkUrl}>Barcode demux benchmark report</a>
              <a href={publicCrisprUrl}>Public CRISPR benchmark report</a>
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
        </div>
      </section>

      <section className="section decision-section" aria-label="DotMatch adoption guide">
        <div className="section-heading">
          <h2>Use it when assignment choices matter.</h2>
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
            Use the source install until the public package channels finish
            publication. Current distribution: source install and repository
            release artifacts, with Bioconda review tracked in PR #65367.
            Channel availability is claimed only after the package appears on
            that channel.
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
            <div className="link-stack">
              <a href={repoUrl}>Open GitHub</a>
              <a href={packagingUrl}>Packaging notes</a>
              <a href={biocondaPrUrl}>Bioconda recipe PR</a>
            </div>
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

      <section className="section scope-section">
        <div className="section-heading">
          <h2>What is validated, early, or out of scope.</h2>
          <p>
            CRISPR guide counting is the strongest public evidence today. Other
            surfaces are useful, but the site keeps smoke tests and future
            distribution work separate from the primary evidence.
          </p>
        </div>
        <div className="scope-layout">
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

          <div className="scope-side">
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
            <div className="scope-notes">
              {evidenceNotes.map(([label, value, detail]) => (
                <article key={label}>
                  <span>{label}</span>
                  <strong>{value}</strong>
                  <p>{detail}</p>
                </article>
              ))}
            </div>
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

const repoUrl = "https://github.com/dnncha/dotmatch";
const citationUrl = `${repoUrl}/blob/main/CITATION.cff`;
const methodsUrl = `${repoUrl}/blob/main/docs/methods-and-citation.md`;
const packagingUrl = `${repoUrl}/blob/main/docs/packaging.md`;
const benchmarksUrl = `${repoUrl}/blob/main/docs/benchmarks/README.md`;
const evidenceGalleryUrl = `${repoUrl}/blob/main/docs/evidence-gallery/README.md`;
const publicCrisprUrl = `${repoUrl}/blob/main/docs/benchmarks/public_crispr/README.md`;
const barcodeBenchmarkUrl = `${repoUrl}/blob/main/docs/benchmarks/barcode_demux/README.md`;
const panelDesignUrl = `${repoUrl}/blob/main/docs/barcode-panel-design.md`;
const panelBenchmarkUrl = `${repoUrl}/blob/main/docs/benchmarks/barcode_panel_design/README.md`;
const biocondaUrl = "https://anaconda.org/bioconda/dotmatch";

const proof = [
  ["Guide counts", "CRISPR libraries", "FASTQ reads become guide-by-sample count tables."],
  ["Barcode splits", "fixed-position indexes", "Inline barcode reads can be assigned, split, and reviewed."],
  ["Panel certificates", "design, check, simulate", "Designed barcodes ship with machine-checkable safety files."],
  ["No guessing", "ambiguity stays visible", "Reads that fit multiple targets are reported as ambiguous."]
];

const decisionCards = [
  {
    title: "Use DotMatch when you have",
    items: [
      "fixed-window barcode FASTQs",
      "CRISPR guide-counting reads",
      "barcode panels to design or certify",
      "known primer, panel, or whitelist targets",
      "feature-barcode or amplicon slices"
    ]
  },
  {
    title: "DotMatch gives you",
    items: [
      "one assignment per read",
      "explicit ambiguous and unmatched reads",
      "panel safety certificates and lab exports",
      "unsafe one-edit correction warnings",
      "HTML, TSV, JSON, and workflow-ready files"
    ]
  },
  {
    title: "Do not use DotMatch for",
    items: [
      "genome alignment or variant calling",
      "basecalling or UMI entropy generation",
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
    title: "Panel designers",
    body: "Generate or audit barcode panels with exact DotMatch assignment certificates, plate layouts, and sample-sheet exports."
  },
  {
    title: "CRISPR screen users",
    body: "Count guides from FASTQ/FASTQ.gz into MAGeCK-compatible matrices, with exact, rescued, ambiguous, and unmatched reads in the QC."
  },
  {
    title: "Methods reviewers",
    body: "Reproduce the exact commands, inspect raw CSVs, and compare reports with full target-by-target or Edlib checks."
  }
];

const workflowStatusRows = [
  ["Barcode panel design", "Good fit", "Design, optimize, check, simulate, lay out, and export panels for known-target assignment."],
  ["CRISPR guide counting", "Good fit", "Guide-by-sample counts, QC summaries, and MAGeCK-compatible output."],
  ["Inline barcode demux", "Good fit", "Fixed-position barcodes, split FASTQs, unmatched reads, and ambiguous reads."],
  ["Barcode troubleshooting", "Good fit", "Window scans, barcode-library checks, and top-unmatched summaries."],
  ["Target-library audit", "Good fit", "Duplicate and near-neighbor checks before one-edit correction."],
  ["Classic BCL demux", "Limited", "Use Illumina BCL Convert for production run-folder conversion."],
  ["Genome alignment", "Use another tool", "DotMatch does not produce SAM/BAM/CIGAR or call variants."]
];

const workflowChoiceRows = [
  ["Design or certify a barcode panel", "DotMatch panel"],
  ["Count CRISPR guides from a fixed window", "DotMatch"],
  ["Split fixed-position inline barcodes", "DotMatch"],
  ["Find why a barcode lane is mostly unassigned", "DotMatch barcode troubleshooting"],
  ["Trim general adapters", "Cutadapt-style tools"],
  ["Map reads to a genome or transcriptome", "Bowtie2, BWA, or minimap2-style tools"],
  ["Analyze CRISPR screen phenotypes", "MAGeCK or another downstream analysis tool"]
];

const evidenceNotes = [
  ["Assignment rule", "fast mode matches full scan", "The indexed search is tested against a check of every target for the same settings."],
  ["Input", "known short targets", "Guides, barcodes, primers, panels, and whitelist-style sequences."],
  ["Repository", "command line and Python", "Core code, bindings, tests, reports, file formats, and benchmark tables."]
];

const commands = [
  "dotmatch panel design --n 96 --length 16 --preset illumina-inline-strict --seed 42 --out-dir panel_96x16",
  "dotmatch panel check panel_96x16/barcodes.tsv --k 1 --metric hamming --out-dir panel_check",
  "dotmatch barcode autopsy --barcodes barcodes.tsv --reads pooled.fastq.gz --scan-starts 0:12 --k-values 0,1 --out-dir autopsy",
  "dotmatch barcode infer --barcodes barcodes.tsv --reads pooled.fastq.gz --scan-starts 0:30 --sample-reads 100000 --out inference.tsv",
  "dotmatch barcode demux --barcodes barcodes.tsv --reads pooled.fastq.gz --barcode-start 1 --barcode-length auto --k 1 --metric hamming --max-correction-qual 20 --out-dir demuxed --report report.html",
  "dotmatch crispr-count --library guides.csv --samples samples.tsv --guide-start 23 --guide-length 19 --k 1 --metric levenshtein --indel-window 1 --out counts.mageck.tsv --summary qc.json",
  "dotmatch assay run assay.toml",
  "dotmatch validate --targets guides.tsv --reads sample.fastq.gz --target-start 23 --target-length 19 --k 1 --indel-window 1 --oracle edlib --sample 100000"
];

const autopsyArtifacts = [
  ["report.html", "HTML summary to open first"],
  ["findings.tsv", "likely offset, rescue, and collision issues"],
  ["offset_scan.tsv", "candidate barcode windows ranked by assignment rate"],
  ["correction_safety.tsv", "whether one-edit rescue can mix barcodes"],
  ["top_unmatched.tsv", "high-count unassigned barcode sequences"],
  ["provenance.json", "run record: commands, versions, thresholds, and output files"]
];

const autopsyFindings = [
  ["wrong offset", "Detects a likely leading base, primer scar, or shifted barcode window."],
  ["unsafe correction", "Shows barcode pairs or clusters that make one-mismatch rescue unsafe."],
  ["ambiguous collision", "Keeps reads that match multiple barcodes out of forced assignments."],
  ["unmatched classes", "Separates low-complexity, distant, reverse-complement, and quality-gated failures."]
];

const panelOutputs = [
  ["barcodes.tsv", "auditable barcode table, not just sequence IDs"],
  ["panel_summary.json", "machine-checkable safety certificate"],
  ["ambiguous_error_spheres.tsv", "queries that would create ambiguity"],
  ["target_safety.tsv", "per-barcode nearest-neighbor and risk status"],
  ["plate_layout.tsv", "96-well or 384-well operational layout"],
  ["SampleSheet.csv", "lab-ready sample-sheet template"]
];

const panelChecks = [
  ["Checked correction rules", "All possible barcode variants are checked up to k=2; larger edit distances are refused."],
  ["Sequence filters", "GC, homopolymer, repeats, forbidden motifs, ambiguous bases, and reverse-complement traps."],
  ["Context checks", "Optional flanks expose cross-boundary homopolymers, motifs, and boundary risks."],
  ["Simulation", "Simple error models estimate unique, ambiguous, none, invalid, and false assignment rates."]
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
const assignmentWorkflowImage = `${basePath}/dotmatch-read-assignment.svg`;
const panelCertificateImage = `${basePath}/dotmatch-panel-certificate.png`;

export default function Home() {
  return (
    <main>
      <header className="site-header">
        <a className="brand" href="#top" aria-label="DotMatch home">
          <span className="brand-mark" />
          DotMatch
        </a>
        <nav aria-label="Primary navigation">
          <a href="#barcode-qc">Barcode QC</a>
          <a href="#panel-design">Panel design</a>
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
            Design panels. Count guides. Split barcodes. See what failed.
          </p>
          <p className="hero-text">
            DotMatch works when the expected short DNA sequences are already
            known: CRISPR guides, inline barcodes, primers, panels, feature
            tags, or whitelist entries. It designs barcode panels, writes count
            or split outputs, and keeps ambiguous, unmatched, and invalid reads
            visible.
          </p>
          <p className="hero-note">
            <strong>Use it after FASTQs exist.</strong>{" "}
            DotMatch does not replace BCL Convert, basecallers, genome aligners,
            or general adapter trimming. It is for known short-DNA target
            assignment.
          </p>
          <div className="hero-actions">
            <a href="#barcode-qc" className="button primary">
              Troubleshoot Barcodes
            </a>
            <a href="#panel-design" className="button secondary">
              Design Panels
            </a>
            <a href={benchmarksUrl} className="button secondary">
              Read Examples
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
            <span>v0.1.2</span>
            <span>known-target assignment</span>
          </div>
          <figure className="hero-visual">
            <img
              src={assignmentWorkflowImage}
              alt="DotMatch workflow showing FASTQ reads and a target table, a fixed read slice, assignment outcomes, and output files"
              decoding="async"
              fetchPriority="high"
            />
            <figcaption>
              FASTQ reads become unique, ambiguous, none, and invalid outcomes,
              with QC tables and reports kept beside panel, count, or split
              outputs.
            </figcaption>
          </figure>
          <div className="metric-grid">
            <div>
              <strong>4</strong>
              <span>assignment outcomes, including ambiguous and invalid reads</span>
            </div>
            <div>
              <strong>8</strong>
              <span>barcode troubleshooting and panel-safety checks in the examples</span>
            </div>
            <div>
              <strong>1.37M</strong>
              <span>reads/s on the checked exact-prefix barcode example</span>
            </div>
            <div>
              <strong>0</strong>
              <span>forced assignments for reads DotMatch reports as ambiguous</span>
            </div>
          </div>
          <div className="sequence-rail" aria-hidden="true">
            {Array.from({ length: 64 }).map((_, i) => (
              <span key={i} className={i % 7 === 0 ? "hot" : i % 5 === 0 ? "cool" : ""} />
            ))}
          </div>
        </div>
      </section>

      <section className="evidence-strip" aria-label="DotMatch workflow summary">
        {proof.map(([label, value, detail]) => (
          <article key={label}>
            <strong>{label}</strong>
            <span>{value}</span>
            <p>{detail}</p>
          </article>
        ))}
      </section>

      <section id="panel-design" className="section panel-design-section">
        <div className="section-heading">
          <h2>Design barcode panels with the assignment rules attached.</h2>
          <p>
            DotMatch panel design creates barcode sets, checks them under the
            same assignment rules used later, and writes safety files a pipeline
            can inspect. It does not hide ambiguous rescue, and it refuses edit
            distances it cannot check exactly.
          </p>
        </div>
        <div className="panel-design-layout">
          <figure className="panel-design-visual">
            <img
              src={panelCertificateImage}
              alt="A lab bench scene with an abstract panel safety report, 96-well plate, and barcode strips"
              decoding="async"
            />
          </figure>
          <article className="panel-command">
            <span className="card-label">Panel lifecycle</span>
            <pre><code>{`dotmatch panel design \\
  --n 96 \\
  --length 16 \\
  --preset illumina-inline-strict \\
  --min-hamming-distance 5 \\
  --min-levenshtein-distance 4 \\
  --seed 42 \\
  --out-dir panel_96x16

dotmatch panel check panel_96x16/barcodes.tsv \\
  --k 1 \\
  --metric hamming \\
  --out-dir panel_check`}</code></pre>
            <p>
              The certificate preserves DotMatch outcomes: unique, ambiguous,
              none, and invalid. DotMatch currently checks all possible barcode
              variants through k=2.
            </p>
            <div className="link-stack compact">
              <a href={panelDesignUrl}>Read panel design docs</a>
              <a href={panelBenchmarkUrl}>Open checked panel example</a>
            </div>
          </article>
        </div>
        <div className="panel-output-grid" aria-label="Panel design outputs">
          {panelOutputs.map(([name, detail]) => (
            <article key={name}>
              <code>{name}</code>
              <p>{detail}</p>
            </article>
          ))}
        </div>
        <div className="panel-check-grid" aria-label="Panel safety checks">
          {panelChecks.map(([name, detail]) => (
            <article key={name}>
              <span>{name}</span>
              <p>{detail}</p>
            </article>
          ))}
        </div>
      </section>

      <section id="barcode-qc" className="section autopsy-section">
        <div className="section-heading">
          <h2>Find the barcode window before you trust the split.</h2>
          <p>
            Barcode demultiplexing can fail because the window is shifted, the
            barcode list has near-neighbors, or one-mismatch rescue would mix
            samples. DotMatch scans the window, checks the barcode table, and
            shows the reads that were not assigned.
          </p>
        </div>
        <div className="autopsy-layout">
          <article className="autopsy-command">
            <span className="card-label">Common command</span>
            <pre><code>{`dotmatch barcode autopsy \\
  --barcodes barcodes.tsv \\
  --reads pooled.fastq.gz \\
  --scan-starts 0:12 \\
  --k-values 0,1 \\
  --out-dir autopsy`}</code></pre>
            <p>
              The command writes the report, window scan, barcode safety table,
              top-unmatched table, and run record into one directory.
            </p>
          </article>
          <div className="artifact-grid" aria-label="Barcode QC outputs">
            {autopsyArtifacts.map(([name, detail]) => (
              <article key={name}>
                <code>{name}</code>
                <p>{detail}</p>
              </article>
            ))}
          </div>
        </div>
        <div className="finding-list" aria-label="Barcode diagnosis examples">
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
          <h2>Reports a lab can actually read.</h2>
          <p>
            Open the HTML report first. Keep the TSV and JSON files for the
            pipeline, notebook, MultiQC page, or methods supplement.
          </p>
        </div>
        <div className="report-preview" aria-label="DotMatch report outcome preview">
          <div className="report-copy">
            <h3>Every read keeps its assignment reason.</h3>
            <p>
              DotMatch separates ambiguous rescue, wrong windows, invalid
              slices, and true no-match reads. That makes it easier to decide
              whether the assay spec is wrong, the barcode list is unsafe, or
              the sample needs to be rerun.
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
          <h2>Examples with commands and raw tables.</h2>
          <p>
            The repository includes fixed-window barcode and CRISPR examples
            with data sources, commands, comparator notes, and raw result files.
            On the repeated public Yusa CRISPR rows, DotMatch Hamming k=1
            processed about 331k reads/s using about 28.7 MB peak memory.
          </p>
        </div>
        <div className="benchmark-grid">
          <article className="benchmark-card">
            <div className="chart-copy">
              <span className="card-label">Evidence gallery</span>
              <h3>See what clean and suspicious runs look like.</h3>
              <p>
                The gallery links public benchmark reports, barcode autopsy
                HTML, findings tables, raw files, and exact commands for
                known-good lanes and diagnostic failure patterns.
              </p>
            </div>
            <div className="link-stack compact">
              <a href={evidenceGalleryUrl}>Open evidence gallery</a>
              <a href={`${repoUrl}/blob/main/docs/evidence-gallery/report-zoo/README.md`}>Open report examples</a>
            </div>
          </article>

          <article className="benchmark-card">
            <div className="chart-copy">
              <span className="card-label">Public CRISPR example</span>
              <h3>CRISPR guide counting.</h3>
              <p>
                Five 100k-record/sample repeats compare DotMatch, MAGeCK, and
                guide-counter on the same public guide-counting example. Exact,
                Hamming, and Levenshtein settings are reported separately.
                Edlib validation checks 2,000 reads with zero mismatches.
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
              <h3>One-edit matching without scanning every guide.</h3>
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
              <h3>Comparator counts are shown beside DotMatch output.</h3>
              <p>
                MAGeCK and guide-counter are useful references for familiar
                CRISPR workflows. DotMatch also checks assignment behavior
                against a full target-by-target scan and Edlib.
              </p>
            </div>
            <AgreementChart rows={agreementRows} />
          </article>
        </div>
      </section>

      <section className="section decision-section" aria-label="DotMatch use guide">
        <div className="section-heading">
          <h2>Use it when assignment choices matter.</h2>
          <p>
            Most DotMatch jobs start as FASTQ reads and a target table. The
            point is not only speed; it is making corrected, ambiguous, and
            unmatched reads visible enough to review.
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
            A typical CRISPR guide-counting run takes reads and a guide library
            and writes a guide-by-sample count matrix plus QC files.
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
          <h2>Install from Bioconda.</h2>
          <p>
            DotMatch 0.1.2 is available from Bioconda for Linux and Intel macOS
            Conda environments. Source builds remain available for other
            platforms.
          </p>
        </div>
        <div className="launch-grid">
          <article className="launch-card">
            <span className="card-label">Bioconda</span>
            <h3>Create an environment and run one command.</h3>
            <pre><code>{`mamba create -n dotmatch -c conda-forge -c bioconda dotmatch=0.1.2
conda activate dotmatch
dotmatch dist ACGT AGGT`}</code></pre>
            <div className="link-stack">
              <a href={biocondaUrl}>Open Bioconda package</a>
              <a href={repoUrl}>Build from source</a>
              <a href={packagingUrl}>Packaging notes</a>
            </div>
          </article>

          <article id="cite" className="launch-card">
            <span className="card-label">Cite it</span>
            <h3>Use the release citation and a matching methods sentence.</h3>
            <p>
              If DotMatch helps an analysis, cite the software release. The
              methods note has short text for CRISPR guide counting,
              one-edit Levenshtein rescue, and Hamming-only comparisons.
            </p>
            <div className="link-stack">
              <a href={citationUrl}>CITATION.cff</a>
              <a href={methodsUrl}>Methods and citation notes</a>
            </div>
          </article>

          <article className="launch-card">
            <span className="card-label">Check the data</span>
            <h3>Read the example before quoting numbers.</h3>
            <p>
              The public CRISPR benchmark is a Yusa-style guide-counting
              example with checked-in rows and assignment validation. Broader
              comparisons need their own datasets and commands.
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
          <h2>Who uses it.</h2>
          <p>
            DotMatch is for people who need short reads assigned to a known
            target list and want the uncertain reads kept visible.
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
          <h2>Where DotMatch fits.</h2>
          <p>
            DotMatch is a focused assignment tool. It is useful when the target
            sequences and read window are known; it is not a genome aligner,
            variant caller, or BCL Convert replacement.
          </p>
        </div>
        <div className="scope-layout">
          <div className="status-table" role="table" aria-label="DotMatch workflow maturity">
            <div role="row" className="table-head">
              <span>Workflow</span>
              <span>Fit</span>
              <span>What it gives you</span>
            </div>
            {workflowStatusRows.map(([workflow, status, evidence]) => (
              <div role="row" key={workflow}>
                <span data-label="Workflow">{workflow}</span>
                <span data-label="Fit">{status}</span>
                <span data-label="What it gives you">{evidence}</span>
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
            DotMatch is a small C/Python tool with a CLI and Python bindings.
            Runs can write count matrices, FASTQ splits, QC tables, assignment
            diagnostics, library checks, validation summaries, and static HTML
            reports.
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
        <h2>Use it when the uncertain reads matter.</h2>
        <p>
          DotMatch is built for fixed-window FASTQ assignment where ambiguous,
          unmatched, and invalid reads should be visible alongside the counts.
        </p>
        <a className="button primary" href="#benchmarks">
          Read Examples
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

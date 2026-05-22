const repoUrl = "https://github.com/dnncha/dotmatch";
const citationUrl = `${repoUrl}/blob/main/CITATION.cff`;
const methodsUrl = `${repoUrl}/blob/main/docs/methods-and-citation.md`;
const packagingUrl = `${repoUrl}/blob/main/docs/packaging.md`;
const benchmarksUrl = `${repoUrl}/blob/main/docs/benchmarks/README.md`;
const evidenceGalleryUrl = `${repoUrl}/blob/main/docs/evidence-gallery/README.md`;
const scientificClaimsUrl = `${repoUrl}/blob/main/docs/scientific-claims.md`;
const publicCrisprUrl = `${repoUrl}/blob/main/docs/benchmarks/public_crispr/README.md`;
const ampliconPanelUrl = `${repoUrl}/blob/main/docs/benchmarks/amplicon_panel/README.md`;
const perturbSeqUrl = `${repoUrl}/blob/main/docs/benchmarks/perturb_seq/README.md`;
const barcodeBenchmarkUrl = `${repoUrl}/blob/main/docs/benchmarks/barcode_demux/README.md`;
const panelDesignUrl = `${repoUrl}/blob/main/docs/barcode-panel-design.md`;
const panelBenchmarkUrl = `${repoUrl}/blob/main/docs/benchmarks/barcode_panel_design/README.md`;
const nextflowExampleUrl = `${repoUrl}/tree/main/examples/workflows/nextflow`;
const nfcoreExampleUrl = `${repoUrl}/tree/main/examples/workflows/nf-core`;
const biocondaPrUrl = "https://github.com/bioconda/bioconda-recipes/pull/65367";
const biocondaUrl = "https://anaconda.org/bioconda/dotmatch";

const proof = [
  ["Guide counts", "screen reads", "FASTQ to guide-by-sample counts."],
  ["Barcode splits", "inline barcodes", "Split reads, then inspect what did not split."],
  ["Panel checks", "barcode panels", "Keep the safety check with the panel."],
  ["No guessing", "ambiguous stays ambiguous", "A read that fits two targets is not forced into one."]
];

const decisionCards = [
  {
    title: "Good for",
    items: [
      "fixed-window barcode FASTQs",
      "CRISPR guide counts",
      "barcode panel design",
      "primer or whitelist checks",
      "feature-barcode windows"
    ]
  },
  {
    title: "You get",
    items: [
      "one assignment per read",
      "ambiguous reads kept separate",
      "HTML reports",
      "TSV and JSON outputs",
      "warnings for unsafe rescue"
    ]
  },
  {
    title: "Not for",
    items: [
      "genome alignment",
      "variant calling",
      "basecalling",
      "cell/UMI processing",
      "CRISPR hit calling"
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
    body: "Find the barcode window, split reads, and see why reads were left out."
  },
  {
    title: "Sequencing cores",
    body: "Check shifted barcode windows, barcode collisions, and unsafe one-mismatch rescue."
  },
  {
    title: "Panel designers",
    body: "Design barcode panels and keep the assignment checks with the panel files."
  },
  {
    title: "CRISPR screen users",
    body: "Count guides from FASTQ into MAGeCK-style tables, with QC beside the counts."
  },
  {
    title: "Methods reviewers",
    body: "Check the commands, raw tables, and validation notes."
  }
];

const workflowStatusRows = [
  ["Barcode panel design", "Good fit", "Design, check, simulate, and export barcode panels."],
  ["CRISPR guide counting", "Good fit", "Guide-by-sample counts and MAGeCK-style output."],
  ["Inline barcode demux", "Good fit", "Split FASTQs and report unmatched or ambiguous reads."],
  ["Barcode troubleshooting", "Good fit", "Scan windows and show likely failure modes."],
  ["Target-library audit", "Good fit", "Find duplicates and near-neighbors before rescue."],
  ["Classic BCL demux", "Limited", "Use Illumina BCL Convert for production run-folder conversion."],
  ["Genome alignment", "Use another tool", "DotMatch does not write SAM/BAM/CIGAR or call variants."]
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
  ["Assignment rule", "index matches scan", "Checked against exhaustive scan for the same settings."],
  ["Input", "known short targets", "Guides, barcodes, primers, panels, or whitelists."],
  ["Repository", "C, CLI, Python", "Code, tests, reports, schemas, and benchmark tables."]
];

const commands = [
  {
    surface: "Python workflow install",
    command: "dotmatch panel design --n 96 --length 16 --preset illumina-inline-strict --seed 42 --out-dir panel_96x16"
  },
  {
    surface: "Python workflow install",
    command: "dotmatch panel check panel_96x16/barcodes.tsv --k 1 --metric hamming --out-dir panel_check"
  },
  {
    surface: "Python workflow install",
    command: "dotmatch barcode autopsy --barcodes barcodes.tsv --reads pooled.fastq.gz --scan-starts 0:12 --k-values 0,1 --out-dir autopsy"
  },
  {
    surface: "Python workflow install",
    command: "dotmatch barcode infer --barcodes barcodes.tsv --reads pooled.fastq.gz --scan-starts 0:30 --sample-reads 100000 --out inference.tsv"
  },
  {
    surface: "Python workflow install",
    command: "dotmatch barcode demux --barcodes barcodes.tsv --reads pooled.fastq.gz --barcode-start 1 --barcode-length auto --k 1 --metric hamming --max-correction-qual 20 --out-dir demuxed --report report.html"
  },
  {
    surface: "Native CLI or Python install",
    command: "dotmatch crispr-count --library guides.csv --samples samples.tsv --guide-start 23 --guide-length 19 --k 1 --metric levenshtein --indel-window 1 --out counts.mageck.tsv --summary qc.json"
  },
  {
    surface: "Python workflow install",
    command: "dotmatch assay run assay.toml"
  },
  {
    surface: "Native CLI or Python install",
    command: "dotmatch validate --targets guides.tsv --reads sample.fastq.gz --target-start 23 --target-length 19 --k 1 --indel-window 1 --oracle edlib --sample 100000"
  }
];

const autopsyArtifacts = [
  ["report.html", "open this first"],
  ["findings.tsv", "likely offset, rescue, and collision issues"],
  ["offset_scan.tsv", "candidate barcode windows ranked by assignment rate"],
  ["correction_safety.tsv", "whether one-edit rescue can mix barcodes"],
  ["top_unmatched.tsv", "high-count unassigned barcode sequences"],
  ["provenance.json", "commands and versions"]
];

const autopsyFindings = [
  ["wrong offset", "The barcode window may be shifted."],
  ["unsafe correction", "One-mismatch rescue may mix samples."],
  ["ambiguous collision", "A read fits more than one barcode."],
  ["unmatched classes", "Common no-match patterns are listed separately."]
];

const panelOutputs = [
  ["barcodes.tsv", "barcode table"],
  ["panel_summary.json", "safety summary"],
  ["ambiguous_error_spheres.tsv", "ambiguous rescue examples"],
  ["target_safety.tsv", "nearest-neighbor checks"],
  ["plate_layout.tsv", "plate layout"],
  ["SampleSheet.csv", "sample-sheet template"]
];

const panelChecks = [
  ["Exact check", "Error spheres are checked up to k=2."],
  ["Sequence filters", "GC, homopolymers, repeats, motifs, and N bases."],
  ["Context checks", "Optional flanks catch boundary problems."],
  ["Simulation", "Estimate unique, ambiguous, none, invalid, and false calls."]
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
const nfcoreArticImage = `${basePath}/dotmatch-nfcore-artic-flow.png`;
const nextflowCrisprImage = `${basePath}/dotmatch-nextflow-crispr-flow.png`;
const guideCaptureImage = `${basePath}/dotmatch-10x-guide-capture-flow.png`;

const realWorkflowExamples = [
  {
    title: "ARTIC V3 primer check in nf-core viralrecon",
    label: "nf-core example",
    image: nfcoreArticImage,
    alt: "Workflow diagram showing R1 FASTQ and ARTIC V3 primers flowing through DotMatch primer-start assignment into QC outputs",
    question: "Do the R1 reads start with the expected ARTIC primer?",
    body: "A public viralrecon R1 file is checked against 80 ARTIC V3 primers. DotMatch uses bases 1-22 of R1 and reports k=0 and k=1 matches.",
    command: `dotmatch count \\
  --targets examples/amplicon_panel/data/artic_v3_primers_len22.tsv \\
  --reads examples/amplicon_panel/data/nfcore_viralrecon_sample1_R1.subsample20000.fastq.gz \\
  --target-start 0 \\
  --target-length 22 \\
  --k 1 \\
  --metric hamming`,
    outputs: ["counts.tsv", "assignments.tsv", "summary.json", "sample_qc.tsv"],
    boundary:
      "Primer-start QC only. Consensus, primer trimming, variant calling, and clinical calls happen elsewhere.",
    links: [
      ["Amplicon evidence", ampliconPanelUrl],
      ["nf-core-style module notes", nfcoreExampleUrl]
    ]
  },
  {
    title: "Nextflow / nf-core-style CRISPR guide counting",
    label: "Nextflow example",
    image: nextflowCrisprImage,
    alt: "Workflow diagram showing a Nextflow DSL2 DotMatch CRISPR guide-counting process with sample FASTQs and guide library inputs",
    question: "Where does guide counting go in a Nextflow screen?",
    body: "The example process stages samples.tsv, FASTQs, and a guide library. DotMatch writes MAGeCK-style counts plus sample QC.",
    command: `nextflow run examples/workflows/nextflow/main.nf \\
  -c examples/workflows/nextflow/nextflow.config`,
    outputs: ["counts.mageck.tsv", "summary.json", "sample_qc.tsv", "assay_report.html"],
    boundary:
      "This is a local module example. MAGeCK hit calling and biology interpretation come after.",
    links: [
      ["Nextflow example", nextflowExampleUrl],
      ["Public CRISPR benchmark", publicCrisprUrl]
    ]
  },
  {
    title: "10x CRISPR guide-capture fixed-window assignment",
    label: "public guide-capture lane",
    image: guideCaptureImage,
    alt: "Workflow diagram showing a 10x CRISPR Guide Capture R2 read window assigned by DotMatch to guide targets with QC outputs",
    question: "Can the guide window be counted before single-cell analysis?",
    body: "A public 10x Guide Capture R2 file is checked at start 63, length 19. DotMatch reports guide counts and per-read assignments.",
    command: `dotmatch count \\
  --targets examples/perturb_seq/data/crispr_guides.tsv \\
  --reads examples/perturb_seq/data/1k_CRISPR_5p_gemx_crispr_S1_L001_R2.subsample20000.fastq.gz \\
  --target-start 63 \\
  --target-length 19 \\
  --k 1 \\
  --metric hamming`,
    outputs: ["counts.tsv", "assignments.tsv", "summary.json", "sample_qc.tsv"],
    boundary:
      "Per-read guide assignment only. Cell barcodes, UMIs, expression matrices, and perturbation calls stay elsewhere.",
    links: [
      ["Guide-capture evidence", perturbSeqUrl],
      ["Evidence gallery", evidenceGalleryUrl]
    ]
  }
];

export default function Home() {
  return (
    <main>
      <header className="site-header">
        <a className="brand" href="#top" aria-label="DotMatch home">
          <span className="brand-mark" />
          DotMatch
        </a>
        <nav aria-label="Primary navigation">
          <a href="#real-workflows">Examples</a>
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
            Use DotMatch when you already know the short sequences you expect:
            guides, barcodes, primers, panels, feature tags, or whitelists. It
            counts or splits reads and shows what did not fit.
          </p>
          <p className="hero-note">
            <strong>After FASTQ.</strong>{" "}
            Not a basecaller, aligner, BCL Convert replacement, or adapter
            trimmer.
          </p>
          <div className="hero-actions">
            <a href="#real-workflows" className="button primary">
              See Real Workflows
            </a>
            <a href="#barcode-qc" className="button secondary">
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
            <span>v0.1.3</span>
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
              Each read is unique, ambiguous, unmatched, or invalid. The report
              stays next to the count or split files.
            </figcaption>
          </figure>
          <div className="metric-grid">
            <div>
              <strong>4</strong>
              <span>outcomes: unique, ambiguous, unmatched, invalid</span>
            </div>
            <div>
              <strong>8</strong>
              <span>barcode and panel checks in the examples</span>
            </div>
            <div>
              <strong>1.37M</strong>
              <span>reads/s in the exact-prefix barcode example</span>
            </div>
            <div>
              <strong>0</strong>
              <span>forced calls for ambiguous reads</span>
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

      <section id="real-workflows" className="section real-workflows-section">
        <div className="section-heading">
          <h2>Real workflow examples.</h2>
          <p>
            Three examples from the repo. Each starts with FASTQ and a known
            target list. DotMatch handles the fixed read window and writes files
            a pipeline can pick up.
          </p>
        </div>
        <div className="real-workflow-list">
          {realWorkflowExamples.map((example) => (
            <article key={example.title} className="real-workflow-card">
              <figure className="real-workflow-visual">
                <img src={example.image} alt={example.alt} decoding="async" loading="eager" />
              </figure>
              <div className="real-workflow-copy">
                <span className="card-label">{example.label}</span>
                <h3>{example.title}</h3>
                <p className="workflow-question">{example.question}</p>
                <p>{example.body}</p>
                <pre><code>{example.command}</code></pre>
                <div className="output-chip-row" aria-label={`${example.title} outputs`}>
                  {example.outputs.map((output) => (
                    <code key={output}>{output}</code>
                  ))}
                </div>
                <p className="boundary-note">{example.boundary}</p>
                <div className="link-stack compact">
                  {example.links.map(([label, href]) => (
                    <a key={label} href={href}>
                      {label}
                    </a>
                  ))}
                </div>
              </div>
            </article>
          ))}
        </div>
      </section>

      <section id="panel-design" className="section panel-design-section">
        <div className="section-heading">
          <h2>Design barcodes, then check the rescue rules.</h2>
          <p>
            A barcode panel is only useful if rescue is safe. DotMatch designs
            panels, checks nearest neighbors, and writes files a lab or pipeline
            can keep.
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
            <span className="card-label">Python workflow install</span>
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
              The check records unique, ambiguous, unmatched, and invalid
              outcomes. Exact error-sphere checks are supported through k=2.
              Use a PyPI or source Python install for these panel commands.
            </p>
            <div className="link-stack compact">
              <a href={panelDesignUrl}>Read panel design docs</a>
              <a href={panelBenchmarkUrl}>Open panel design gate</a>
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
            If a split looks wrong, start with the window. DotMatch scans likely
            starts, checks the barcode list, and shows the reads it could not
            assign.
          </p>
        </div>
        <div className="autopsy-layout">
          <article className="autopsy-command">
            <span className="card-label">Python workflow install</span>
            <pre><code>{`dotmatch barcode autopsy \\
  --barcodes barcodes.tsv \\
  --reads pooled.fastq.gz \\
  --scan-starts 0:12 \\
  --k-values 0,1 \\
  --out-dir autopsy`}</code></pre>
            <p>
              One directory: report, window scan, barcode safety, top unmatched,
              and provenance. Use a PyPI or source Python install for barcode
              troubleshooting commands.
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
          <h2>A report you can open during QC.</h2>
          <p>
            Open the HTML first. Keep the TSV and JSON files for the pipeline,
            notebook, MultiQC page, or methods section.
          </p>
        </div>
        <div className="report-preview" aria-label="DotMatch report outcome preview">
          <div className="report-copy">
            <h3>Every read keeps a reason.</h3>
            <p>
              Ambiguous rescue, wrong windows, invalid slices, and no-match
              reads are separate. That makes the next decision clearer.
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
          <h2>Examples you can rerun.</h2>
          <p>
            The repo includes public FASTQ examples, commands, raw tables, and
            comparator notes. The Yusa CRISPR repeats are included so the
            numbers can be checked.
          </p>
        </div>
        <div className="benchmark-grid">
          <article className="benchmark-card">
            <div className="chart-copy">
              <span className="card-label">Evidence gallery</span>
              <h3>Clean runs. Suspicious runs. Same format.</h3>
              <p>
                Public reports, autopsy HTML, findings tables, raw artifacts,
                and the commands that made them.
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
                Five repeats compare DotMatch, MAGeCK, and guide-counter on the
                same public guide-counting example. Exact, Hamming, and
                Levenshtein runs are kept separate.
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
                On the Yusa rows, one-edit Levenshtein checks about 2.8
                candidate guides per read, not the whole 87,437-guide library.
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
                In the repeated Yusa runs, DotMatch exact and Hamming lanes sit
                around 28.7 MB peak memory.
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
              <h3>Counts are compared, not hand-waved.</h3>
              <p>
                MAGeCK, guide-counter, exhaustive scan, and Edlib are used where
                they answer the right question.
              </p>
            </div>
            <AgreementChart rows={agreementRows} />
          </article>
        </div>
      </section>

      <section className="section decision-section" aria-label="DotMatch use guide">
        <div className="section-heading">
          <h2>Use it for short known targets.</h2>
          <p>
            Most jobs are just FASTQ plus a target table. The important part is
            seeing exact, rescued, ambiguous, and unmatched reads separately.
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
            Reads in. Guide library in. Counts and QC out.
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
              Ambiguous reads are not counted into a guide or barcode. They stay
              visible.
            </p>
          </article>
        </div>
      </section>

      <section id="install" className="section launch-section">
        <div className="section-heading">
          <h2>Install the native CLI from Bioconda.</h2>
          <p>
            Bioconda provides the native `dotmatch` command plus C
            header/library artifacts on published platforms. Use a PyPI or
            source Python install for `dotmatch assay`, `dotmatch barcode`, and
            `dotmatch panel` workflow commands.
          </p>
        </div>
        <div className="launch-grid">
          <article className="launch-card">
            <span className="card-label">Native CLI</span>
            <h3>Create a Bioconda env.</h3>
            <pre><code>{`conda create -n dotmatch -c conda-forge -c bioconda dotmatch=0.1.2
conda activate dotmatch
dotmatch --version
dotmatch dist ACGT AGGT`}</code></pre>
            <div className="link-stack">
              <a href={biocondaUrl}>Open Bioconda package</a>
              <a href={repoUrl}>Open GitHub</a>
              <a href={packagingUrl}>Packaging notes</a>
              <a href={biocondaPrUrl}>Bioconda recipe PR</a>
            </div>
          </article>

          <article className="launch-card">
            <span className="card-label">Python workflow layer</span>
            <h3>Use PyPI or source installs.</h3>
            <p>
              The workflow namespaces for assay specs, barcode autopsy, panel
              design, HTML reports, and Workbench-backed runs are not part of
              the native Bioconda package.
            </p>
            <div className="link-stack">
              <a href={packagingUrl}>Read package boundaries</a>
              <a href={repoUrl}>Install from source</a>
            </div>
          </article>

          <article id="cite" className="launch-card">
            <span className="card-label">Cite it</span>
            <h3>Cite the release.</h3>
            <p>
              The methods note has short text for CRISPR counting, Levenshtein
              rescue, and Hamming-only runs.
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
              The public CRISPR benchmark is a Yusa guide-counting example with
              checked rows and validation notes.
            </p>
            <div className="link-stack">
              <a href={publicCrisprUrl}>Public CRISPR benchmark report</a>
              <a href={scientificClaimsUrl}>Scientific claims and boundaries</a>
              <a href="#benchmarks">Review benchmark summary</a>
            </div>
          </article>
        </div>
      </section>

      <section id="use-cases" className="section use-cases">
        <div className="section-heading">
          <h2>Who uses it.</h2>
          <p>
            For people who need short reads assigned to a known list, with the
            uncertain reads left visible.
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
            DotMatch is a small assignment tool. It works when the read window
            and target sequences are known.
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
            Native CLI commands cover fixed-window assignment and validation.
            Python workflow installs add barcode autopsy, panel design, assay
            specs, and HTML reports.
          </p>
        </div>
        <div className="terminal" aria-label="DotMatch commands">
          <div className="terminal-bar">
            <span />
            <span />
            <span />
          </div>
          {commands.map(({ surface, command }) => (
            <code key={command}>
              <span className="terminal-surface">{surface}</span>
              <span>$</span> {command}
            </code>
          ))}
        </div>
      </section>

      <section className="section final-cta">
        <h2>Use it when the uncertain reads matter.</h2>
        <p>
          Fixed-window FASTQ assignment, with ambiguous and unmatched reads kept
          in view.
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

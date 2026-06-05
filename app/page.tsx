const repoUrl = "https://github.com/dnncha/dotmatch";
const doi = "10.5281/zenodo.20541629";
const doiUrl = `https://doi.org/${doi}`;
const citationUrl = `${repoUrl}/blob/main/CITATION.cff`;
const methodsUrl = `${repoUrl}/blob/main/docs/methods-and-citation.md`;
const packagingUrl = `${repoUrl}/blob/main/docs/packaging.md`;
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
const initialBiocondaPrUrl = "https://github.com/bioconda/bioconda-recipes/pull/65367";
const biocondaUrl = "https://anaconda.org/bioconda/dotmatch";

const proof = [
  ["CRISPR screens", "guide counts", "Count guides from FASTQ and write MAGeCK-style tables."],
  ["Barcode demux", "inline barcodes", "Split reads and keep unmatched or ambiguous reads visible."],
  ["Panel design", "barcode sets", "Design panels and check whether error correction is safe."],
  ["QC reports", "HTML + tables", "Share run summaries without losing the raw TSV and JSON files."],
  ["Fast paths", "benchmarked", "Public examples include commands, raw data, and graphs."]
];

const decisionCards = [
  {
    title: "Supported use cases",
    items: [
      "CRISPR guide counting",
      "inline barcode demultiplexing",
      "barcode panel design",
      "primer or whitelist checks",
      "feature-barcode reads"
    ]
  },
  {
    title: "Outputs",
    items: [
      "one assignment per read",
      "ambiguity tables",
      "HTML reports",
      "TSV and JSON outputs",
      "unsafe-rescue diagnostics"
    ]
  },
  {
    title: "Out of scope",
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
  ["Target list", "the guide, barcode, primer, or whitelist sequences you expect"],
  ["One mismatch", "allow one substituted base, with no insertions or deletions"],
  ["GuideCounter mode", "one-mismatch Hamming by default, exact-only with --exact-match"],
  ["One edit", "allow one substitution, insertion, or deletion"],
  ["Ambiguous", "a read matches more than one target and is not forced into either one"],
  ["Peak RSS", "peak memory use"],
  ["Edlib validation", "checked against an independent edit-distance implementation"]
];

const audienceCards = [
  {
    title: "Barcode assay owners",
    body: "Estimate barcode windows, demultiplex reads, and review ambiguous or unmatched classes."
  },
  {
    title: "Sequencing cores",
    body: "Audit shifted barcode windows, barcode collisions, and unsafe one-mismatch rescue."
  },
  {
    title: "Panel designers",
    body: "Design barcode panels, check collision risk, and export lab-ready files."
  },
  {
    title: "CRISPR screen users",
    body: "Count guides from FASTQ into MAGeCK-style tables, with QC beside the counts."
  },
  {
    title: "Methods reviewers",
    body: "Inspect commands, raw tables, comparison notes, and validation checks."
  }
];

const workflowStatusRows = [
  ["Barcode panel design", "Supported", "Design, check, simulate, and export barcode panels."],
  ["CRISPR guide counting", "Supported", "Guide-by-sample counts and MAGeCK-style output."],
  ["Inline barcode demux", "Supported", "Split FASTQs and report unmatched or ambiguous reads."],
  ["Barcode troubleshooting", "Supported", "Scan candidate windows and summarize failure modes."],
  ["Target-library audit", "Supported", "Identify duplicates and near-neighbors before rescue."],
  ["Classic BCL demux", "Limited", "Use Illumina BCL Convert for production run-folder conversion."],
  ["Genome alignment", "Out of scope", "DotMatch does not write SAM/BAM/CIGAR or call variants."]
];

const workflowChoiceRows = [
  ["Design or check a barcode panel", "DotMatch panel"],
  ["Count CRISPR guides from a fixed window", "DotMatch"],
  ["Split fixed-position inline barcodes", "DotMatch"],
  ["Diagnose a low-yield barcode lane", "DotMatch barcode troubleshooting"],
  ["Trim general adapters", "Cutadapt-style tools"],
  ["Map reads to a genome or transcriptome", "Bowtie2, BWA, or minimap2-style tools"],
  ["Analyze CRISPR screen phenotypes", "MAGeCK or another downstream analysis tool"]
];

const evidenceNotes = [
  ["Matching rule", "index matches scan", "Checked against exhaustive scan for the same settings."],
  ["Input", "known short sequences", "Guides, barcodes, primers, panels, or whitelists."],
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
    surface: "Native CLI",
    command: "dotmatch guide-counter count --input sample.fastq.gz --library guides.tsv --output guide_counts"
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
  ["report.html", "HTML summary"],
  ["findings.tsv", "offset, rescue, and collision diagnostics"],
  ["offset_scan.tsv", "candidate barcode windows ranked by assignment rate"],
  ["correction_safety.tsv", "whether one-edit rescue can mix barcodes"],
  ["top_unmatched.tsv", "high-count unassigned barcode sequences"],
  ["provenance.json", "commands and versions"]
];

const autopsyFindings = [
  ["wrong offset", "Candidate barcode windows show a shifted assignment peak."],
  ["unsafe correction", "One-mismatch rescue creates cross-barcode compatibility."],
  ["ambiguous collision", "A read fits more than one barcode."],
  ["unmatched classes", "Common no-match patterns are listed separately."]
];

const panelOutputs = [
  ["barcodes.tsv", "barcode table"],
  ["panel_summary.json", "collision-risk summary"],
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
  ["ambiguous", "compatible with multiple targets", "written to ambiguity output"],
  ["none", "outside the configured edit radius", "sent to unmatched diagnostics"],
  ["invalid", "window could not be extracted", "recorded in QC"]
];

const throughputRows = [
  { label: "DotMatch exact k=0", value: 887206.3, tone: "green" },
  { label: "DotMatch Hamming k=1", value: 754902.5, tone: "green" },
  { label: "guide-counter one mismatch", value: 184061.4, tone: "blue" },
  { label: "MAGeCK exact count", value: 127848.7, tone: "gray" },
  { label: "DotMatch Levenshtein k=1", value: 6635.0, tone: "green" }
] as const;

const memoryRows = [
  { label: "guide-counter one mismatch", value: 528.7, tone: "blue" },
  { label: "MAGeCK exact count", value: 152.9, tone: "gray" },
  { label: "DotMatch exact k=0", value: 113.4, tone: "green" },
  { label: "DotMatch Hamming k=1", value: 49.3, tone: "green" },
  { label: "DotMatch Levenshtein k=1", value: 113.5, tone: "green" }
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

const benchmarkFigures = [
  {
    title: "Public CRISPR throughput",
    label: "CRISPR comparison",
    image: `${basePath}/benchmarks/crispr_comparison_throughput.svg`,
    alt: "CRISPR comparison throughput graph for DotMatch, MAGeCK, and guide-counter public datasets",
    body: "Exact, Hamming, and Levenshtein lanes are separated so the chart compares the same biological question rather than mixing rescue rules.",
    href: `${repoUrl}/blob/main/docs/benchmarks/crispr_comparison/README.md`
  },
  {
    title: "Hamming k2/k3 comparator",
    label: "Bowtie 1 comparison",
    image: `${basePath}/benchmarks/crispr_hamming_k23_comparison.svg`,
    alt: "Hamming k2 and k3 fixed-window CRISPR comparator graph comparing DotMatch with Bowtie 1",
    body: "Hamming k2/k3 rows are separated from GuideCounter claims. Bowtie 1 is used here for same-strand fixed-window Hamming comparisons.",
    href: `${repoUrl}/blob/main/docs/benchmarks/crispr_comparison/README.md`
  },
  {
    title: "Apple Metal GPU lane",
    label: "Experimental GPU",
    image: `${basePath}/benchmarks/gpu_crispr_metal_speedup.svg`,
    alt: "Public CRISPR GPU benchmark graph comparing DotMatch CPU indexed Hamming assignment with the Apple Metal packed Hamming lane",
    body: "The public CRISPR row includes FASTQ parsing, guide-window extraction, GPU dispatch, readback, and count aggregation. It is evidence for productizing GPU work, not a default production claim.",
    href: `${repoUrl}/blob/main/docs/benchmarks/gpu/README.md`
  },
  {
    title: "Synthetic GPU stress test",
    label: "GPU comparison",
    image: `${basePath}/benchmarks/gpu_metal_speedup.svg`,
    alt: "Synthetic GPU benchmark graph comparing Apple Metal brute-force Hamming assignment against DotMatch CPU indexed assignment",
    body: "The synthetic lane shows where a GPU can help when the target set is large and every read-target pair is packable A/C/G/T Hamming k=1.",
    href: `${repoUrl}/blob/main/docs/benchmarks/gpu/README.md`
  },
  {
    title: "Inline barcode demux",
    label: "Barcode comparison",
    image: `${basePath}/benchmarks/barcode_demux_throughput.svg`,
    alt: "Barcode demultiplexing throughput graph for DotMatch, Cutadapt, and exact hash splitter rows",
    body: "The barcode graph uses a public exact-prefix SRP009896 lane. It compares fixed-position barcode matching, anchored Cutadapt demux, and a simple exact-prefix baseline.",
    href: barcodeBenchmarkUrl
  },
  {
    title: "Barcode memory",
    label: "Resource use",
    image: `${basePath}/benchmarks/barcode_demux_peak_memory.svg`,
    alt: "Peak memory graph for the public barcode demultiplexing benchmark",
    body: "Memory is shown beside throughput because a demux tool is only useful if it keeps routine lanes practical on shared machines.",
    href: barcodeBenchmarkUrl
  },
  {
    title: "Repeated CRISPR throughput",
    label: "Public repeatability",
    image: `${basePath}/benchmarks/public_crispr_repeated_throughput.svg`,
    alt: "Repeated public CRISPR benchmark graph showing throughput by tool and edit-distance lane",
    body: "Repeated rows make run-to-run variation visible before anyone quotes a single throughput number.",
    href: publicCrisprUrl
  },
  {
    title: "Repeated CRISPR memory",
    label: "Resource use",
    image: `${basePath}/benchmarks/public_crispr_repeated_peak_memory.svg`,
    alt: "Repeated public CRISPR benchmark graph showing peak memory by tool and edit-distance lane",
    body: "Memory stays beside speed so the comparison is useful for real shared workstations and CI machines.",
    href: publicCrisprUrl
  },
  {
    title: "Verified work per read",
    label: "Algorithm shape",
    image: `${basePath}/benchmarks/public_crispr_repeated_verified_candidates.svg`,
    alt: "Public CRISPR repeated benchmark graph showing verified candidate guides per read",
    body: "This is the reason the indexed path matters: one-edit Levenshtein checks a small candidate set instead of scanning the whole guide library.",
    href: publicCrisprUrl
  }
] as const;

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
      "Workflow integration only: DotMatch emits the count matrix and QC; MAGeCK, BAGEL, drugZ, or CERES handle screen-level inference.",
    links: [
      ["Nextflow example", nextflowExampleUrl],
      ["Public CRISPR benchmark", publicCrisprUrl]
    ]
  },
  {
    title: "10x CRISPR guide-capture fixed-window check",
    label: "single-guide public lane",
    image: guideCaptureImage,
    alt: "Workflow diagram showing a 10x CRISPR Guide Capture R2 read window assigned by DotMatch to guide targets with QC outputs",
    question: "Can the guide sequence window be checked before single-cell analysis?",
    body: "A public 10x Guide Capture R2 file is checked at start 63, length 19. The example is a fixed-window single-guide extraction and count check.",
    command: `dotmatch count \\
  --targets examples/perturb_seq/data/crispr_guides.tsv \\
  --reads examples/perturb_seq/data/1k_CRISPR_5p_gemx_crispr_S1_L001_R2.subsample20000.fastq.gz \\
  --target-start 63 \\
  --target-length 19 \\
  --k 1 \\
  --metric hamming`,
    outputs: ["counts.tsv", "assignments.tsv", "summary.json", "sample_qc.tsv"],
    boundary:
      "Single-guide fixed-window check only. Cell barcodes, UMIs, expression matrices, multi-guide calls, and perturbation calls stay elsewhere.",
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
          <a href="#real-workflows">Workflows</a>
          <a href="#barcode-qc">Barcode QC</a>
          <a href="#panel-design">Panel design</a>
          <a href="#benchmarks">Benchmarks</a>
          <a href="#use-cases">Use cases</a>
          <a href="#install">Install</a>
          <a href="#cite">Cite</a>
          <a href={repoUrl}>GitHub</a>
        </nav>
        <a className="header-cta" href={repoUrl}>Source</a>
        <div className="mobile-header-actions" aria-label="Quick navigation">
          <a href="#real-workflows">Examples</a>
          <a href="#install">Install</a>
        </div>
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
              See workflows
            </a>
            <a href={repoUrl} className="button secondary">
              GitHub
            </a>
          </div>
          <div className="hero-link-row" aria-label="Secondary DotMatch links">
            <a href="#barcode-qc">Barcode QC</a>
            <a href="#panel-design">Panel Design</a>
            <a href="#benchmarks">Benchmarks</a>
            <a href="#install">Install</a>
          </div>
        </div>
        <div className="hero-panel" aria-label="DotMatch benchmark summary">
          <div className="panel-topline">
            <span>v0.1.7 archived DOI</span>
            <span>FASTQ sequence matching</span>
          </div>
          <figure className="hero-visual">
            <img
              src={assignmentWorkflowImage}
              alt="DotMatch workflow showing FASTQ reads and a target table, a fixed read slice, assignment outcomes, and output files"
              decoding="async"
              fetchPriority="high"
            />
            <figcaption>
              Every read gets a clear outcome. The report stays next to the count
              or demultiplexed files.
            </figcaption>
          </figure>
          <div className="metric-grid">
            <div>
              <strong>4</strong>
              <span>read outcomes tracked in every run</span>
            </div>
            <div>
              <strong>8</strong>
              <span>barcode and panel checks in the examples</span>
            </div>
            <div>
              <strong>2.4M</strong>
              <span>reads/s in the exact-prefix barcode example</span>
            </div>
            <div>
              <strong>0</strong>
              <span>forced calls when two targets fit</span>
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
          <h2>Workflows you can actually run.</h2>
          <p>
            Each example starts with FASTQ and a target list. DotMatch handles the
            read window, writes useful files, and leaves a report behind for QC.
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
          <h2>Design barcode panels with correction in mind.</h2>
          <p>
            A barcode set should stay reliable after sequencing errors. DotMatch
            designs panels, checks nearest neighbors, and exports files a lab or
            pipeline can use.
          </p>
        </div>
        <div className="panel-design-layout">
          <figure className="panel-design-visual">
            <img
              src={panelCertificateImage}
              alt="A lab bench scene with an abstract barcode panel collision report, 96-well plate, and barcode strips"
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
              The check shows whether one- or two-error correction could mix
              samples. Use a source Python install, or PyPI after the tagged
              release is visible, for these panel commands.
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
        <div className="panel-check-grid" aria-label="Panel collision checks">
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
          <h2>Find the barcode window before demultiplexing.</h2>
          <p>
            DotMatch scans likely barcode positions, checks the barcode list, and
            shows unmatched or ambiguous reads before you split the data.
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
              One directory contains the report, window scan, correction checks,
              top unmatched sequences, and provenance. Use a source Python
              install, or PyPI after the tagged release is visible, for barcode
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
          <h2>A QC report people can read.</h2>
          <p>
            Open the HTML report during review. Keep the TSV and JSON files for
            the pipeline, notebook, MultiQC page, or methods section.
          </p>
        </div>
        <div className="report-preview" aria-label="DotMatch report outcome preview">
          <div className="report-copy">
            <h3>No silent guessing.</h3>
            <p>
              Ambiguous matches, wrong windows, invalid slices, and no-match reads
              are separated so the next decision is obvious.
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
          <h2>Benchmarks with the details included.</h2>
          <p>
            The repo includes public FASTQ examples, commands, raw tables,
            generated graphs, and notes on what each tool was asked to do. CPU
            indexed assignment remains the production baseline. GPU rows are an
            experimental evidence lane and must match CPU output before speed
            matters.
          </p>
        </div>
        <div className="benchmark-reader-guide" aria-label="How to read DotMatch benchmarks">
          <article>
            <span>Compare like with like</span>
            <p>Exact matching, one-mismatch matching, deeper Hamming rescue, and one-edit rescue are reported separately.</p>
          </article>
          <article>
            <span>Check the comparison</span>
            <p>MAGeCK, guide-counter, Bowtie 1, Cutadapt, Edlib, and hash baselines are used only where their behavior matches the question.</p>
          </article>
          <article>
            <span>Treat GPU as experimental</span>
            <p>Metal rows must match CPU outputs before speed is considered useful.</p>
          </article>
        </div>
        <div className="benchmark-figure-grid">
          {benchmarkFigures.map((figure) => (
            <figure key={figure.title} className="benchmark-figure">
              <a href={figure.href} aria-label={`Open source report for ${figure.title}`}>
                <img src={figure.image} alt={figure.alt} loading="lazy" decoding="async" />
              </a>
              <figcaption>
                <span className="card-label">{figure.label}</span>
                <strong>{figure.title}</strong>
                <p>{figure.body}</p>
              </figcaption>
            </figure>
          ))}
        </div>
        <div className="benchmark-grid">
          <article className="benchmark-card">
            <div className="chart-copy">
              <span className="card-label">Evidence gallery</span>
              <h3>Clean runs. Suspicious runs. Same format.</h3>
              <p>
                Public reports, troubleshooting HTML, findings tables, raw artifacts,
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
              <span className="card-label">GuideCounter-compatible lane</span>
              <h3>CRISPR guide counting with familiar output files.</h3>
              <p>
                Repeated public rows compare DotMatch exact against MAGeCK
                exact, DotMatch Hamming k=1 against guide-counter one-mismatch,
                and DotMatch Levenshtein k=1 as its own indel-rescue lane. The
                compatibility command writes <code>counts.txt</code>,{" "}
                <code>extended-counts.txt</code>, and <code>stats.txt</code>{" "}
                for teams that already have GuideCounter-shaped downstream steps.
              </p>
            </div>
            <div className="link-stack compact">
              <a href={publicCrisprUrl}>Public CRISPR benchmark report</a>
              <a href={`${repoUrl}/blob/main/docs/benchmarks/crispr_comparison/README.md`}>CRISPR comparison report</a>
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
              <h3>One-edit matching without brute force.</h3>
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
              <h3>CRISPR guide counting keeps memory low.</h3>
              <p>
                In the repeated Yusa 100k rows, DotMatch Hamming sits around
                49 MB peak memory; exact and Levenshtein sit around 113 MB.
              </p>
            </div>
            <HorizontalBarChart
              rows={memoryRows}
              unit="MB"
              axisLabel="Peak memory, lower is better"
              scale="linear"
            />
          </article>

          <article className="benchmark-card">
            <div className="chart-copy">
              <span className="card-label">Count agreement</span>
              <h3>Count agreement is checked against the right references.</h3>
              <p>
                MAGeCK, guide-counter, exhaustive scan, and Edlib are used where
                their behavior matches the reported comparison.
              </p>
            </div>
            <AgreementChart rows={agreementRows} />
          </article>
        </div>
      </section>

      <section className="section decision-section" aria-label="DotMatch use guide">
        <div className="section-heading">
          <h2>When DotMatch is the right tool.</h2>
          <p>
            Use DotMatch when each read should be compared with a known list of
            guides, barcodes, primers, features, or whitelist sequences.
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
        <div className="translation-grid" aria-label="Plain-language glossary for DotMatch terms">
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
          <h2>One CRISPR run, from FASTQ to counts.</h2>
          <p>
            Reads in. Guide library in. Counts and QC out. Use
            <code>guide-counter count</code> compatibility when existing
            scripts expect GuideCounter-style files.
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
          <article className="example-card">
            <span className="card-label">GuideCounter-compatible</span>
            <pre><code>{`dotmatch guide-counter count \\
  --input plasmid.fastq.gz treatment.fastq.gz \\
  --samples plasmid treatment \\
  --library guides.tsv \\
  --output guide_counts`}</code></pre>
            <div className="output-list" aria-label="GuideCounter-compatible outputs">
              <code>guide_counts.counts.txt</code>
              <span>guide x sample count matrix</span>
              <code>guide_counts.extended-counts.txt</code>
              <span>guide_type plus counts</span>
              <code>guide_counts.stats.txt</code>
              <span>mapped fraction and guide-class means</span>
            </div>
            <p className="boundary-note">
              Defaults match the one-mismatch/no-indel Hamming lane with
              automatic multi-offset detection. Add <code>--exact-match</code>{" "}
              for exact-only counting.
            </p>
          </article>
          <article className="ambiguity-example">
            <span className="card-label">Why ambiguity matters</span>
            <pre><code>{`Read:    ACGTACGT
Guide A: ACGTACGA   distance 1
Guide B: ACGTACGC   distance 1

Some tools may pick or double-count.
	DotMatch reports: ambiguous`}</code></pre>
            <p>
              Ambiguous reads are not counted into a guide or barcode. They stay
              visible for review.
            </p>
          </article>
        </div>
      </section>

      <section id="install" className="section launch-section">
        <div className="section-heading">
          <h2>Install the stable package or build the current release.</h2>
          <p>
            Bioconda currently publishes DotMatch 0.1.4. Newer features in this
            branch should be installed from source until the matching package
            version passes public channel smoke tests.
          </p>
        </div>
        <div className="launch-grid">
          <article className="launch-card">
            <span className="card-label">Published Bioconda package</span>
            <h3>Create a verified 0.1.4 env.</h3>
            <pre><code>{`conda create -n dotmatch -c conda-forge -c bioconda dotmatch=0.1.4
conda activate dotmatch
dotmatch --version
dotmatch dist ACGT AGGT`}</code></pre>
            <div className="link-stack">
              <a href={biocondaUrl}>Open Bioconda package</a>
              <a href={repoUrl}>Open GitHub</a>
              <a href={packagingUrl}>Packaging notes</a>
              <a href={initialBiocondaPrUrl}>Initial Bioconda PR</a>
            </div>
          </article>

          <article className="launch-card">
            <span className="card-label">Release candidate</span>
            <h3>Use source installs for the newest workflows.</h3>
            <p>
              GuideCounter-compatible mode, new evidence reports, and release
              candidate packaging require a source install until the tagged
              release is verified across PyPI, Bioconda, containers, and DOI
              metadata.
            </p>
            <div className="link-stack">
              <a href={packagingUrl}>Read package boundaries</a>
              <a href={repoUrl}>Install from source</a>
            </div>
          </article>

          <article id="cite" className="launch-card">
            <span className="card-label">Cite it</span>
            <h3>Cite the archived release.</h3>
            <p>
              Zenodo DOI: <a href={doiUrl}>{doi}</a>. Use the citation file for
              software metadata and the methods note for manuscript wording.
            </p>
            <div className="link-stack">
              <a href={doiUrl}>Open Zenodo DOI</a>
              <a href={citationUrl}>CITATION.cff</a>
              <a href={methodsUrl}>Methods and citation notes</a>
            </div>
          </article>

          <article className="launch-card">
            <span className="card-label">Benchmark evidence</span>
            <h3>Cite performance with the benchmark context.</h3>
            <p>
              The public CRISPR benchmark uses a Yusa guide-counting example with
              checked rows, commands, and validation notes.
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
          <h2>Primary users.</h2>
          <p>
            For groups assigning short read windows to known guide, barcode,
            primer, feature, whitelist, or panel targets while retaining
            ambiguous and unmatched read classes.
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
            DotMatch is a focused FASTQ matching tool. It works when the read
            window and expected sequences are known.
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
            Native CLI commands cover fixed-window matching and validation.
            Python workflow installs add barcode troubleshooting, panel design,
            assay specs, and HTML reports.
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
        <h2>Fast FASTQ matching with honest QC.</h2>
        <p>
          Count guides, split barcodes, design panels, and keep ambiguous,
          unmatched, and invalid reads visible in the output.
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

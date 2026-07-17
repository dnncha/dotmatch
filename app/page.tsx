const repoUrl = "https://github.com/dnncha/dotmatch";
const docsUrl = "https://dotmatch.readthedocs.io/en/latest/";
const pypiUrl = "https://pypi.org/project/dotmatch/";
const biocondaUrl = "https://anaconda.org/bioconda/dotmatch";

const basePath = process.env.NEXT_PUBLIC_BASE_PATH ?? "";
const assignmentWorkflowImage = `${basePath}/dotmatch-read-assignment.svg`;

const structuredData = {
  "@context": "https://schema.org",
  "@graph": [
    {
      "@type": "WebSite",
      "@id": "https://dnncha.github.io/dotmatch/#website",
      name: "DotMatch",
      url: "https://dnncha.github.io/dotmatch",
      description:
        "DotMatch assigns short FASTQ read windows to known DNA targets for CRISPR, barcode, and other fixed-target sequencing workflows."
    },
    {
      "@type": "SoftwareApplication",
      "@id": "https://dnncha.github.io/dotmatch/#software",
      name: "DotMatch",
      applicationCategory: "Bioinformatics software",
      operatingSystem: "Linux, macOS",
      softwareHelp: docsUrl,
      codeRepository: repoUrl,
      license: `${repoUrl}/blob/main/LICENSE`,
      programmingLanguage: ["C", "Python", "R"],
      description:
        "DotMatch compares fixed read windows with known DNA targets and reports unique, ambiguous, unmatched, and invalid reads."
    }
  ]
};

const outcomes = [
  ["unique", "One target is compatible. The read can be counted or written to that target's output."],
  ["ambiguous", "Several targets are compatible. The read is kept out of the unique counts."],
  ["none", "No target is within the selected distance. The read remains available for review."],
  ["invalid", "The requested window could not be extracted. The failure appears in QC."]
] as const;

const failureModes = [
  {
    title: "A read fits more than one target",
    body: "DotMatch reports the ambiguity instead of choosing a target arbitrarily."
  },
  {
    title: "The target window is in the wrong place",
    body: "Offset scans and invalid-window counts help distinguish a bad start position from genuinely unmatched reads."
  },
  {
    title: "Mismatch correction is unsafe",
    body: "Target-library audits show duplicate and neighbouring sequences before corrected reads are counted."
  },
  {
    title: "Unmatched reads repeat",
    body: "Frequent unmatched windows are written to a table so adapter, assay, and sample-sheet problems can be inspected."
  }
] as const;

const workflows = [
  {
    title: "CRISPR guide counting",
    body: "Count known guide windows, write MAGeCK-compatible tables, and retain ambiguous and unmatched reads in QC outputs.",
    link: "Run the CRISPR tutorial",
    href: `${docsUrl}tutorials/crispr-count-first-run.html`
  },
  {
    title: "Inline barcode demultiplexing",
    body: "Split FASTQ reads by fixed-position sample barcodes and inspect offsets, near neighbours, and recurring unmatched sequences.",
    link: "Start with a demultiplexing run",
    href: `${docsUrl}getting-started.html#demultiplex-inline-barcodes`
  },
  {
    title: "Feature tags and fixed targets",
    body: "Assign feature barcodes, guide-capture reads, primers, adapters, amplicon starts, or another finite target list.",
    link: "Read the command reference",
    href: `${docsUrl}command-reference.html`
  },
  {
    title: "Barcode panel checks",
    body: "Design a panel or check whether the selected correction radius can create ambiguous or incorrect assignments.",
    link: "Open the panel guide",
    href: `${docsUrl}barcode-panel-design.html`
  }
] as const;

const documentationLinks = [
  ["Getting started", `${docsUrl}getting-started.html`],
  ["Command reference", `${docsUrl}command-reference.html`],
  ["Output schemas", `${docsUrl}schemas.html`],
  ["Python API", `${docsUrl}streaming-api.html`],
  ["Benchmarks", `${docsUrl}benchmarks/README.html`],
  ["Methods and citation", `${docsUrl}methods-and-citation.html`]
] as const;

export default function Home() {
  return (
    <>
      <script
        type="application/ld+json"
        dangerouslySetInnerHTML={{ __html: JSON.stringify(structuredData) }}
      />
      <header className="site-header">
        <a className="brand" href="#top" aria-label="DotMatch home">
          <span className="brand-mark" aria-hidden="true" />
          <span>DotMatch</span>
        </a>
        <nav aria-label="Primary navigation">
          <a href="#how-it-works">How it works</a>
          <a href="#workflows">Workflows</a>
          <a href="#documentation">Documentation</a>
          <a href="#install">Install</a>
          <a href={repoUrl}>GitHub</a>
        </nav>
      </header>

      <main id="top">
        <section className="hero" aria-labelledby="hero-title">
          <div className="hero-copy">
            <p className="positioning">Known-target assignment from FASTQ</p>
            <h1 id="hero-title">Keep ambiguous reads out of your counts.</h1>
            <p className="hero-lede">
              DotMatch compares a fixed read window with the guides, barcodes,
              feature tags, primers, or other short DNA targets you expect.
            </p>
            <p className="hero-text">
              Each read is recorded as unique, ambiguous, unmatched, or invalid.
              The result is a set of ordinary FASTQ, TSV, JSON, and HTML files that
              can be reviewed by a person or passed into a pipeline.
            </p>
            <div className="hero-actions">
              <a className="button primary" href="#install">Install DotMatch</a>
              <a className="button secondary" href={`${docsUrl}getting-started.html`}>Read the guide</a>
            </div>
          </div>
          <div className="hero-panel" aria-label="DotMatch read outcomes">
            <figure className="assignment-figure">
              <img
                src={assignmentWorkflowImage}
                alt="FASTQ reads and a known target table are compared at a fixed read window and written to counts, split FASTQs, QC tables, and reports"
                decoding="async"
                fetchPriority="high"
              />
            </figure>
            <div className="outcome-grid" aria-label="Per-read assignment outcomes">
              {outcomes.map(([label, detail]) => (
                <article key={label} className={`outcome-card ${label}`}>
                  <strong>{label}</strong>
                  <p>{detail}</p>
                </article>
              ))}
            </div>
          </div>
        </section>

        <section id="how-it-works" className="section failure-section" aria-labelledby="how-title">
          <div className="section-heading">
            <p className="section-kicker">What stays visible</p>
            <h2 id="how-title">The difficult reads are part of the result.</h2>
            <p>
              A known target list does not make every assignment safe. DotMatch
              writes the cases that need attention instead of hiding them inside a
              count matrix.
            </p>
          </div>
          <div className="failure-grid">
            {failureModes.map((mode) => (
              <article key={mode.title}>
                <h3>{mode.title}</h3>
                <p>{mode.body}</p>
              </article>
            ))}
          </div>
        </section>

        <section id="workflows" className="section audience-section" aria-labelledby="workflow-title">
          <div className="section-heading">
            <p className="section-kicker">Workflows</p>
            <h2 id="workflow-title">One assignment rule, several assay types.</h2>
            <p>
              DotMatch is deliberately narrow: it works when the expected short
              sequences are known and the read window can be stated explicitly.
            </p>
          </div>
          <div className="audience-grid">
            {workflows.map((workflow) => (
              <article key={workflow.title}>
                <h3>{workflow.title}</h3>
                <p>{workflow.body}</p>
                <a href={workflow.href}>{workflow.link}</a>
              </article>
            ))}
          </div>
        </section>

        <section id="example" className="section evidence-section" aria-labelledby="example-title">
          <div className="section-heading">
            <p className="section-kicker">A direct count</p>
            <h2 id="example-title">State the target window and keep the outputs.</h2>
            <p>
              This example assigns a 20-base guide beginning at read position 23.
              Hamming distance allows substitutions but not insertions or deletions.
            </p>
          </div>
          <div className="terminal" aria-label="DotMatch count example">
            <div className="terminal-bar" aria-hidden="true"><span /><span /><span /></div>
            <pre><code>{`dotmatch count \\
  --targets guides.tsv \\
  --reads sample_R1.fastq.gz \\
  --sample-label sample_1 \\
  --target-start 23 \\
  --target-length 20 \\
  --k 1 --metric hamming \\
  --out counts.tsv \\
  --sample-qc sample_qc.tsv \\
  --summary summary.json`}</code></pre>
            <div className="terminal-links">
              <a href={`${docsUrl}getting-started.html`}>Explain this example</a>
              <a href={`${docsUrl}command-reference.html`}>See every command</a>
            </div>
          </div>
        </section>

        <section id="documentation" className="section workflow-section" aria-labelledby="docs-title">
          <div className="section-heading">
            <p className="section-kicker">Documentation</p>
            <h2 id="docs-title">Start with the job you need to do.</h2>
            <p>
              The user guide begins with runnable commands. File formats, APIs,
              benchmarks, and citation details live in separate reference pages.
            </p>
          </div>
          <div className="evidence-links" aria-label="DotMatch documentation links">
            {documentationLinks.map(([label, href]) => (
              <a key={label} href={href}>{label}</a>
            ))}
          </div>
        </section>

        <section id="install" className="section install-section" aria-labelledby="install-title">
          <div className="install-copy">
            <p className="section-kicker">Install locally</p>
            <h2 id="install-title">Run DotMatch where the FASTQ files live.</h2>
            <p>
              Wheels are published for supported Linux and macOS platforms. Conda
              users can use the Bioconda package when the required release is present.
            </p>
          </div>
          <div className="terminal" aria-label="Install DotMatch">
            <div className="terminal-bar" aria-hidden="true"><span /><span /><span /></div>
            <pre><code>{`python3 -m pip install dotmatch
dotmatch --version`}</code></pre>
            <div className="terminal-links">
              <a href={pypiUrl}>PyPI</a>
              <a href={biocondaUrl}>Bioconda</a>
              <a href={`${docsUrl}packaging.html`}>Other installation routes</a>
            </div>
          </div>
        </section>
      </main>

      <footer className="site-footer">
        <span>DotMatch</span>
        <nav aria-label="Footer navigation">
          <a href={docsUrl}>Documentation</a>
          <a href={repoUrl}>GitHub</a>
          <a href={`${docsUrl}methods-and-citation.html`}>Citation</a>
          <a href={`${repoUrl}/blob/main/LICENSE`}>Apache-2.0</a>
        </nav>
      </footer>
    </>
  );
}

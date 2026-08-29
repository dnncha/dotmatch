import { MobileNavigation } from "./mobile-navigation";

const repoUrl = "https://github.com/dnncha/dotmatch";
const docsUrl = "https://dotmatch.readthedocs.io/";
const gettingStartedUrl = "https://dotmatch.readthedocs.io/en/latest/getting-started.html";
const scientificClaimsUrl = `${repoUrl}/blob/main/docs/scientific-claims.md`;
const evidenceGalleryUrl = `${repoUrl}/blob/main/docs/evidence-gallery/README.md`;
const methodsUrl = `${repoUrl}/blob/main/docs/methods-and-citation.md`;
const packagingUrl = `${repoUrl}/blob/main/docs/packaging.md`;
const evaluationUrl = `${repoUrl}/blob/main/docs/bioinformatics-evaluation.md`;
const reviewPacketUrl = `${repoUrl}/blob/main/docs/external-review-packet.md`;
const integrationKitUrl = `${repoUrl}/blob/main/docs/workflow-integration-kit.md`;
const workflowSubmissionsUrl = `${repoUrl}/blob/main/docs/workflow-submissions.md`;
const workflowAdoptionUrl = `${repoUrl}/blob/main/docs/workflow-adoption.json`;
const distributionUrl = `${repoUrl}/blob/main/docs/distribution-release.json`;
const pypiUrl = "https://pypi.org/project/dotmatch/";
const biocondaUrl = "https://anaconda.org/bioconda/dotmatch";
const doiUrl = "https://doi.org/10.5281/zenodo.20541628";

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
        "DotMatch is a deterministic known-target sequencing assignment toolkit for CRISPR guides, inline barcodes, feature tags, primers, and panel targets."
    },
    {
      "@type": "SoftwareApplication",
      "@id": "https://dnncha.github.io/dotmatch/#software",
      name: "DotMatch",
      applicationCategory: "Bioinformatics software",
      operatingSystem: "Linux, macOS",
      softwareHelp: "https://dotmatch.readthedocs.io/",
      codeRepository: repoUrl,
      license: `${repoUrl}/blob/main/LICENSE`,
      programmingLanguage: ["C", "Python", "R"],
      description:
        "DotMatch assigns fixed read windows to known short DNA targets and reports unique, ambiguous, none, and invalid outcomes for auditable sequencing workflows."
    }
  ]
};

const outcomes = [
  ["unique", "Exactly one target is compatible, so the read can be counted or written to the matching output."],
  ["ambiguous", "More than one target is compatible, so DotMatch keeps the read out of forced calls."],
  ["none", "No target is close enough, so the read remains available for unmatched-read review."],
  ["invalid", "The requested read window cannot be extracted, so the failure is visible in QC."]
] as const;

const failureModes = [
  {
    title: "More than one target fits.",
    body:
      "The read is reported as ambiguous instead of being forced into one count."
  },
  {
    title: "The read window is missing.",
    body:
      "Short reads and invalid extraction windows are reported explicitly instead of disappearing."
  },
  {
    title: "Correction could mix samples.",
    body:
      "Target-library checks show when an error-correction setting could assign one sequence to multiple targets."
  },
  {
    title: "No expected target fits.",
    body:
      "Unmatched-read tables preserve recurring sequences for assay, adapter, and off-target review."
  }
] as const;

const workflowSteps = [
  {
    step: "1",
    title: "Provide the sequences you expect.",
    body:
      "Start with a table of guides, inline barcodes, feature tags, primers, panel targets, or whitelist sequences."
  },
  {
    step: "2",
    title: "Choose where to look in each read.",
    body:
      "DotMatch extracts the same configured position from every read and compares it with the expected sequences."
  },
  {
    step: "3",
    title: "Review every outcome.",
    body:
      "Each read is unique, ambiguous, unmatched, or invalid. Counts, FASTQs, QC tables, and reports preserve that decision."
  }
] as const;

const audienceRoutes = [
  {
    title: "Core facilities",
    body:
      "Use DotMatch when sample barcodes, guide libraries, or panel targets need visible ambiguity and unmatched-read review before a result leaves the core.",
    link: "Start with barcode troubleshooting",
    href: `${repoUrl}/blob/main/docs/crispr-qc.md`
  },
  {
    title: "CRISPR screen teams",
    body:
      "Count known guide windows, keep MAGeCK-compatible outputs, and preserve assignment failures for methods review and downstream screen analysis.",
    link: "Run the CRISPR tutorial",
    href: `${repoUrl}/blob/main/docs/tutorials/crispr-count-first-run.md`
  },
  {
    title: "Workflow maintainers",
    body:
      "Wrap stable TSV, JSON, FASTQ, and HTML artifacts in nf-core, Galaxy, Snakemake, MultiQC, or institutional pipeline templates.",
    link: "Use the submission pack",
    href: workflowSubmissionsUrl
  },
  {
    title: "Assay developers",
    body:
      "Design and audit barcode panels, test correction radius safety, and export lab-ready panel records before sequencing starts.",
    link: "Review panel design",
    href: `${repoUrl}/blob/main/docs/barcode-panel-design.md`
  }
] as const;

const contexts = [
  "CRISPR guides",
  "inline barcodes",
  "feature tags",
  "primers / panels",
  "whitelists"
] as const;

const evidenceLinks = [
  ["Bioinformatics evaluation packet", evaluationUrl],
  ["External review packet", reviewPacketUrl],
  ["Validated scope", scientificClaimsUrl],
  ["Evidence gallery", evidenceGalleryUrl],
  ["Methods and citation", methodsUrl],
  ["Packaging notes", packagingUrl]
] as const;

const evaluationItems = [
  {
    title: "Published release",
    body:
      "DotMatch v0.2.2 is available from PyPI and Bioconda. The package pages list the current release and installation files."
  },
  {
    title: "Defined scope",
    body:
      "The documented scope is fixed-position assignment against known short DNA targets. DotMatch is not presented as a genome aligner, basecaller, or downstream analysis package."
  },
  {
    title: "Inspectable outputs",
    body:
      "TSV, JSON, FASTQ, and HTML outputs are documented, alongside methods text, citation metadata, and recorded software versions."
  },
  {
    title: "Reproducible checks",
    body:
      "Correctness, performance, and packaging statements link to the repository artifacts and checks used to support them."
  },
  {
    title: "Workflow examples",
    body:
      "The repository includes nf-core, MultiQC, Galaxy, and Snakemake examples. They are examples maintained here, not claims of accepted upstream integration."
  }
] as const;

const trustFacts = [
  ["Release", "v0.2.2 on PyPI and Bioconda"],
  ["License", "Apache-2.0 open source"],
  ["Data", "Runs locally on your machine"],
  ["Citation", "Archived release with DOI"]
] as const;

const outputTypes = [
  "count tables",
  "split FASTQs",
  "QC tables",
  "JSON summaries",
  "HTML reports"
] as const;

const navigationLinks = [
  ["How it works", "#workflow"],
  ["Reliability", "#failure-modes"],
  ["Uses", "#industry-routes"],
  ["Evidence", "#evidence"],
  ["Install", "#install"]
] as const;

export default function Home() {
  return (
    <>
      <script
        type="application/ld+json"
        dangerouslySetInnerHTML={{ __html: JSON.stringify(structuredData) }}
      />
      <a className="skip-link" href="#main-content">Skip to main content</a>
      <header className="site-header">
        <a className="brand" href="#top" aria-label="DotMatch home">
          <span className="brand-mark" aria-hidden="true" />
          <span>DotMatch</span>
        </a>
        <nav className="desktop-nav" aria-label="Primary navigation">
          {navigationLinks.map(([label, href]) => (
            <a key={href} href={href}>{label}</a>
          ))}
          <a href={docsUrl}>Docs</a>
          <a href={repoUrl}>GitHub</a>
        </nav>
        <MobileNavigation links={navigationLinks} docsUrl={docsUrl} repoUrl={repoUrl} />
      </header>

      <main id="main-content">
        <section id="top" className="hero" aria-labelledby="hero-title">
          <div className="hero-copy">
            <p className="positioning">Open-source known-target read assignment</p>
            <h1 id="hero-title">See which sequencing reads match—and which do not.</h1>
            <p className="hero-lede">
              DotMatch compares a chosen part of each read with a list of short DNA
              sequences you expect. Every read gets one clear outcome: unique,
              ambiguous, no match, or invalid.
            </p>
            <p className="hero-text">
              Use it for CRISPR guides, inline barcodes, feature tags, primers,
              panel targets, and whitelist-style assays. DotMatch is not a genome
              aligner or basecaller.
            </p>
            <div className="hero-actions">
              <a className="button primary" href="#install">Try DotMatch</a>
              <a className="button secondary" href="#workflow">See how it works</a>
            </div>
            <dl className="trust-list" aria-label="DotMatch package facts">
              {trustFacts.map(([term, detail]) => (
                <div key={term}>
                  <dt>{term}</dt>
                  <dd>{detail}</dd>
                </div>
              ))}
            </dl>
          </div>
          <div className="hero-panel" aria-label="DotMatch assignment outcomes">
            <figure className="assignment-figure">
              <img
                src={assignmentWorkflowImage}
                alt="DotMatch workflow showing FASTQ reads, a known target list, a read slice, and unique, ambiguous, none, and invalid outcomes"
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

        <section id="failure-modes" className="section failure-section" aria-labelledby="failure-title">
          <div className="section-heading">
            <p className="section-kicker">Uncertainty stays visible</p>
            <h2 id="failure-title">Uncertain reads stay uncertain.</h2>
            <p>
              A target list does not make every read safe to count. DotMatch keeps
              the difficult cases separate so they can be reviewed instead of hidden.
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

        <section id="workflow" className="section workflow-section" aria-labelledby="workflow-title">
          <div className="section-heading">
            <p className="section-kicker">How it works</p>
            <h2 id="workflow-title">Three steps, one clear result per read.</h2>
            <p>
              The scientific question can vary, but the assignment is simple: provide
              expected sequences, choose a position in the read, and inspect the result.
            </p>
          </div>
          <div className="context-rail" aria-label="Supported known-target assay contexts">
            {contexts.map((context) => (
              <span key={context}>{context}</span>
            ))}
          </div>
          <div className="workflow-grid">
            {workflowSteps.map((item) => (
              <article key={item.step}>
                <span aria-hidden="true">{item.step}</span>
                <h3>{item.title}</h3>
                <p>{item.body}</p>
              </article>
            ))}
          </div>
        </section>

        <section id="industry-routes" className="section audience-section" aria-labelledby="audience-title">
          <div className="section-heading">
            <p className="section-kicker">Where it fits</p>
            <h2 id="audience-title">Built for familiar assay problems.</h2>
            <p>
              DotMatch is useful when the expected short sequences are already known
              and the assignment decision needs to remain inspectable.
            </p>
          </div>
          <div className="audience-grid">
            {audienceRoutes.map((route) => (
              <article key={route.title}>
                <h3>{route.title}</h3>
                <p>{route.body}</p>
                <a href={route.href}>{route.link}</a>
              </article>
            ))}
          </div>
        </section>

        <section id="evidence" className="section evidence-section" aria-labelledby="evidence-title">
          <div className="section-heading">
            <p className="section-kicker">Methods and evidence</p>
            <h2 id="evidence-title">Inspect the methods, limits, and outputs.</h2>
            <p>
              The package documents what it does, what it does not do, how releases
              are checked, and which artifacts support scientific and performance statements.
            </p>
          </div>
          <div className="evidence-layout">
            <article className="evidence-note">
              <h3>Scope, stated plainly</h3>
              <p>
                DotMatch performs deterministic assignment of a fixed read window
                against a known list of short DNA sequences. It keeps ambiguous,
                unmatched, and invalid reads visible alongside unique matches.
              </p>
              <p>
                It is not a general aligner, basecaller, variant caller, UMI pipeline,
                or downstream CRISPR screen-analysis package. The linked scientific
                scope explains the current evidence in detail.
              </p>
            </article>
            <div className="evidence-links" aria-label="Evidence and packaging links">
              {evidenceLinks.map(([label, href]) => (
                <a key={label} href={href}>{label}</a>
              ))}
            </div>
          </div>
        </section>

        <section id="evaluation" className="section evaluation-section" aria-labelledby="evaluation-title">
          <div className="section-heading">
            <p className="section-kicker">Current package status</p>
            <h2 id="evaluation-title">What you can verify today.</h2>
            <p>
              Release records, package channels, output formats, scientific scope,
              and repository checks are available for independent review.
            </p>
          </div>
          <div className="evaluation-layout">
            <ol className="evaluation-list">
              {evaluationItems.map((action) => (
                <li key={action.title}>
                  <strong>{action.title}</strong>
                  <p>{action.body}</p>
                </li>
              ))}
            </ol>
            <aside className="evaluation-links" aria-label="Evaluation and evidence links">
              <a href={evaluationUrl}>Bioinformatics evaluation packet</a>
              <a href={reviewPacketUrl}>External review packet</a>
              <a href={distributionUrl}>Distribution status record</a>
              <a href={workflowAdoptionUrl}>Workflow adoption status</a>
              <a href={docsUrl}>User documentation</a>
              <a href={scientificClaimsUrl}>Scientific scope notes</a>
            </aside>
          </div>
        </section>

        <section id="ecosystem" className="section ecosystem-section" aria-labelledby="ecosystem-title">
          <div className="section-heading">
            <p className="section-kicker">Workflow-friendly outputs</p>
            <h2 id="ecosystem-title">Use ordinary files in ordinary workflows.</h2>
            <p>
              DotMatch writes familiar, inspectable files rather than requiring a
              proprietary project format.
            </p>
          </div>
          <div className="ecosystem-layout">
            <ul className="ecosystem-grid" aria-label="DotMatch output types">
              {outputTypes.map((output) => (
                <li key={output}>{output}</li>
              ))}
            </ul>
            <aside className="ecosystem-note">
              <h3>Pipeline examples</h3>
              <p>
                The repository includes examples for nf-core, MultiQC, Galaxy, and
                Snakemake, plus documented schemas for workflow authors. These are
                maintained DotMatch examples unless an external project records an
                accepted integration.
              </p>
              <a href={workflowSubmissionsUrl}>Workflow submission pack</a>
              <a href={integrationKitUrl}>Workflow integration kit</a>
              <a href={workflowAdoptionUrl}>Recorded integration status</a>
            </aside>
          </div>
        </section>

        <section id="install" className="section install-section" aria-labelledby="install-title">
          <div className="install-copy">
            <p className="section-kicker">Try it locally</p>
            <h2 id="install-title">Start with a two-sequence check.</h2>
            <p>
              Install the command-line and Python package, then run a small distance
              check before working with FASTQ data. Your sequencing files remain on
              your machine.
            </p>
          </div>
          <div className="terminal" aria-label="Install DotMatch">
            <div className="terminal-bar" aria-hidden="true">
              <span />
              <span />
              <span />
            </div>
            <pre><code>{`python -m pip install dotmatch
dotmatch dist ACGT AGGT`}</code></pre>
            <div className="terminal-links">
              <a href={pypiUrl}>PyPI package</a>
              <a href={biocondaUrl}>Bioconda package</a>
              <a href={gettingStartedUrl}>Getting started</a>
              <a href={packagingUrl}>Packaging documentation</a>
              <a href={methodsUrl}>Citation guidance</a>
            </div>
          </div>
        </section>
      </main>

      <footer className="site-footer">
        <span>DotMatch · Apache-2.0</span>
        <nav aria-label="Footer navigation">
          <a href={docsUrl}>Documentation</a>
          <a href={repoUrl}>GitHub</a>
          <a href={scientificClaimsUrl}>Scientific claims</a>
          <a href={methodsUrl}>Methods</a>
          <a href={doiUrl}>DOI</a>
        </nav>
      </footer>
    </>
  );
}

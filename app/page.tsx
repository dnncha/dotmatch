const repoUrl = "https://github.com/dnncha/dotmatch";
const scientificClaimsUrl = `${repoUrl}/blob/main/docs/scientific-claims.md`;
const evidenceGalleryUrl = `${repoUrl}/blob/main/docs/evidence-gallery/README.md`;
const methodsUrl = `${repoUrl}/blob/main/docs/methods-and-citation.md`;
const packagingUrl = `${repoUrl}/blob/main/docs/packaging.md`;
const evaluationUrl = `${repoUrl}/blob/main/docs/bioinformatics-evaluation.md`;
const reviewPacketUrl = `${repoUrl}/blob/main/docs/external-review-packet.md`;
const integrationTargetsUrl = `${repoUrl}/blob/main/docs/integration-targets.json`;
const pilotProgramUrl = `${repoUrl}/blob/main/docs/pilot-program.md`;
const reviewerReadinessUrl = `${repoUrl}/blob/main/docs/reviewer-readiness.json`;
const integrationKitUrl = `${repoUrl}/blob/main/docs/workflow-integration-kit.md`;
const workflowSubmissionsUrl = `${repoUrl}/blob/main/docs/workflow-submissions.md`;
const adoptersUrl = `${repoUrl}/blob/main/docs/adopters/README.md`;
const workflowAdoptionUrl = `${repoUrl}/blob/main/docs/workflow-adoption.json`;
const distributionUrl = `${repoUrl}/blob/main/docs/distribution-release.json`;
const pypiUrl = "https://pypi.org/project/dotmatch/";

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
    title: "A read fits more than one target.",
    body:
      "Ambiguity is reported as ambiguity instead of being hidden inside a count matrix."
  },
  {
    title: "The expected window is wrong.",
    body:
      "Shifted barcode starts, short reads, and invalid extraction windows are surfaced before they become silent losses."
  },
  {
    title: "Correction would mix samples.",
    body:
      "Target-library audits and barcode checks show when rescue settings can create unsafe assignments."
  },
  {
    title: "Unmatched reads carry signal.",
    body:
      "Top-unmatched tables keep recurring off-target, adapter, or assay-design patterns available for review."
  }
] as const;

const workflowSteps = [
  {
    step: "1",
    title: "Declare the known targets.",
    body:
      "Use the guide, inline barcode, feature tag, primer or panel target, or whitelist sequences you expect for the assay."
  },
  {
    step: "2",
    title: "Assign the same read window.",
    body:
      "DotMatch extracts the configured window from each read and evaluates it against the target list under the recorded run settings."
  },
  {
    step: "3",
    title: "Review outcomes with artifacts.",
    body:
      "Counts, split FASTQs, QC tables, unmatched reads, ambiguity rows, summaries, and reports stay connected to the assignment decision."
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
    title: "Package channels",
    body:
      "PyPI and Bioconda are verified for v0.1.8. GHCR and BioContainers have public records, with runtime smoke tests still pending on a host with Docker or another OCI runtime."
  },
  {
    title: "Validated scope",
    body:
      "Public statements are scoped to checked assay lanes, raw artifacts, generated reports, and gate scripts. Broader workflow claims stay out of release copy."
  },
  {
    title: "Output contracts",
    body:
      "TSV, JSON, FASTQ, HTML, methods, citation, and software-version artifacts are documented for workflow systems and reviewer handoff."
  },
  {
    title: "Workflow status",
    body:
      "Local nf-core, MultiQC, Galaxy, and Snakemake examples exist. External workflow integration is recorded only after an accepted public record exists."
  },
  {
    title: "Public use records",
    body:
      "Named labs, projects, organizations, and quotes appear only with approved wording and a public URL."
  }
] as const;

const ecosystemTargets = [
  "nf-core modules",
  "MultiQC module",
  "Galaxy / IUC",
  "Snakemake wrapper",
  "bio.tools record"
] as const;

export default function Home() {
  return (
    <>
      <script
        type="application/ld+json"
        dangerouslySetInnerHTML={{ __html: JSON.stringify(structuredData) }}
      />
      <header className="site-header">
        <a className="brand" href="#top" aria-label="AssayCode home">
          <span className="brand-mark" aria-hidden="true" />
          <span>AssayCode</span>
        </a>
        <nav aria-label="Primary navigation">
          <a href="#failure-modes">Reliability</a>
          <a href="#workflow">Workflow</a>
          <a href="#industry-routes">Use cases</a>
          <a href="#evidence">Evidence</a>
          <a href="#evaluation">Evaluation</a>
          <a href="#ecosystem">Ecosystem</a>
          <a href="#install">Install</a>
          <a href={repoUrl}>GitHub</a>
        </nav>
      </header>

      <main id="top">
        <section className="hero" aria-labelledby="hero-title">
          <div className="hero-copy">
            <p className="positioning">Assay compilation and reliability for known-target sequencing.</p>
            <h1 id="hero-title">Design the assay. Trust the assignment.</h1>
            <p className="hero-lede">
              AssayCode turns known-target assay descriptions into reviewable plans, then uses
              the DotMatch engine to keep every read outcome visible: unique, ambiguous, none, or invalid.
            </p>
            <p className="hero-text">
              Use it when the guide, inline barcode, feature tag, primer or panel
              target, or whitelist sequence is already known and the important question
              is whether each read can be assigned safely.
            </p>
            <div className="hero-actions">
              <a className="button primary" href="#install">Install locally</a>
              <a className="button secondary" href="#evidence">Review evidence</a>
            </div>
          </div>
          <div className="hero-panel" aria-label="AssayCode assignment outcomes">
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
            <p className="section-kicker">Reliability failure modes</p>
            <h2 id="failure-title">The risky reads are the point.</h2>
            <p>
              A known target list does not make every read safe to count. DotMatch is
              built around the cases that need to stay visible.
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
            <p className="section-kicker">One workflow across assay contexts</p>
            <h2 id="workflow-title">The reliability layer is the same.</h2>
            <p>
              CRISPR guides, inline barcodes, feature tags, primers, panel targets,
              and whitelist sequences all reduce to the same auditable assignment
              question.
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
            <p className="section-kicker">Routes into industry workflows</p>
            <h2 id="audience-title">Give each evaluator a next step.</h2>
            <p>
              The people who can spread DotMatch need different entry points:
              core facilities want reliable handoff, screen teams want count
              compatibility, workflow maintainers want stable outputs, and assay
              teams want panel safety.
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
            <p className="section-kicker">Auditable evidence</p>
            <h2 id="evidence-title">Claims stay tied to checks.</h2>
            <p>
              Speed is useful only after the assignment rule is clear. DotMatch keeps
              performance, correctness, and packaging statements scoped to repository
              evidence, public artifacts, and install smoke tests.
            </p>
          </div>
          <div className="evidence-layout">
            <article className="evidence-note">
              <h3>What the homepage can safely claim</h3>
              <p>
                DotMatch is a deterministic known-target assignment system for short
                read windows. It records explicit read outcomes and writes ordinary
                workflow artifacts that can be inspected outside the homepage.
              </p>
              <p>
                Broader claims about alignment, basecalling, downstream screen
                analysis, production demultiplexing replacement, calibrated
                probabilities, or unbounded assay coverage remain outside the public
                validated scope unless the linked evidence says otherwise.
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
            <p className="section-kicker">Bioinformatics evaluation</p>
            <h2 id="evaluation-title">Start from package status and output evidence.</h2>
            <p>
              A serious evaluation should begin with install channels, workflow
              artifacts, claim boundaries, and the exact places where verification
              is still incomplete.
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
              <a href={integrationTargetsUrl}>Integration target tracker</a>
              <a href={scientificClaimsUrl}>Scientific scope notes</a>
            </aside>
          </div>
        </section>

        <section id="ecosystem" className="section ecosystem-section" aria-labelledby="ecosystem-title">
          <div className="section-heading">
            <p className="section-kicker">Workflow ecosystem</p>
            <h2 id="ecosystem-title">Integrations are tracked separately from scientific claims.</h2>
            <p>
              DotMatch has local workflow examples and submission payloads, but
              external integration is recorded only after a public integration is
              accepted or released outside this repository.
            </p>
          </div>
          <div className="ecosystem-layout">
            <ol className="ecosystem-grid" aria-label="Priority workflow integration targets">
              {ecosystemTargets.map((target) => (
                <li key={target}>{target}</li>
              ))}
            </ol>
            <aside className="ecosystem-note">
              <h3>Current rule</h3>
              <p>
                Bioconda package availability is distribution evidence, not an
                accepted workflow integration. Accepted integrations belong in
                the machine-readable workflow status record.
              </p>
              <a href={workflowSubmissionsUrl}>Workflow submission pack</a>
              <a href={integrationKitUrl}>Workflow integration kit</a>
              <a href={pilotProgramUrl}>DotMatch evaluation protocol</a>
              <a href={reviewerReadinessUrl}>Reviewer readiness record</a>
              <a href={methodsUrl}>Methods and citation text</a>
              <a href={adoptersUrl}>Public use record policy</a>
            </aside>
          </div>
        </section>

        <section id="install" className="section install-section" aria-labelledby="install-title">
          <div className="install-copy">
            <p className="section-kicker">Local installation</p>
            <h2 id="install-title">Run assignment reliability checks locally.</h2>
            <p>
              DotMatch is a local command-line and Python package. Keep sequencing
              data on your machine, install the package, and cite the software release
              through the repository citation metadata.
            </p>
          </div>
          <div className="terminal" aria-label="Install DotMatch">
            <div className="terminal-bar" aria-hidden="true">
              <span />
              <span />
              <span />
            </div>
            <pre><code>{`pip install dotmatch
dotmatch --help`}</code></pre>
            <div className="terminal-links">
              <a href={pypiUrl}>PyPI package</a>
              <a href={packagingUrl}>Packaging documentation</a>
              <a href={methodsUrl}>Citation guidance</a>
            </div>
          </div>
        </section>
      </main>

      <footer className="site-footer">
        <span>AssayCode</span>
        <nav aria-label="Footer navigation">
          <a href={repoUrl}>GitHub</a>
          <a href={scientificClaimsUrl}>Scientific claims</a>
          <a href={methodsUrl}>Methods</a>
          <a href={packagingUrl}>Packaging</a>
        </nav>
      </footer>
    </>
  );
}

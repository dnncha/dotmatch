const repoUrl = "https://github.com/dnncha/dotmatch";
const scientificClaimsUrl = `${repoUrl}/blob/main/docs/scientific-claims.md`;
const evidenceGalleryUrl = `${repoUrl}/blob/main/docs/evidence-gallery/README.md`;
const methodsUrl = `${repoUrl}/blob/main/docs/methods-and-citation.md`;
const packagingUrl = `${repoUrl}/blob/main/docs/packaging.md`;
const pypiUrl = "https://pypi.org/project/dotmatch/";

const basePath = process.env.NEXT_PUBLIC_BASE_PATH ?? "";
const assignmentWorkflowImage = `${basePath}/dotmatch-read-assignment.svg`;

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

const contexts = [
  "CRISPR guides",
  "inline barcodes",
  "feature tags",
  "primers / panels",
  "whitelists"
] as const;

const evidenceLinks = [
  ["Claim boundaries", scientificClaimsUrl],
  ["Evidence gallery", evidenceGalleryUrl],
  ["Methods and citation", methodsUrl],
  ["Packaging notes", packagingUrl]
] as const;

export default function Home() {
  return (
    <>
      <header className="site-header">
        <a className="brand" href="#top" aria-label="DotMatch home">
          <span className="brand-mark" aria-hidden="true" />
          <span>DotMatch</span>
        </a>
        <nav aria-label="Primary navigation">
          <a href="#failure-modes">Reliability</a>
          <a href="#workflow">Workflow</a>
          <a href="#evidence">Evidence</a>
          <a href="#install">Install</a>
          <a href={repoUrl}>GitHub</a>
        </nav>
      </header>

      <main id="top">
        <section className="hero" aria-labelledby="hero-title">
          <div className="hero-copy">
            <p className="positioning">Assignment reliability for known-target sequencing assays.</p>
            <h1 id="hero-title">Know which read assignments you can trust.</h1>
            <p className="hero-lede">
              DotMatch assigns fixed read windows to known short DNA targets and keeps
              the outcome visible for every read: unique, ambiguous, none, or invalid.
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
                claim boundary unless the linked evidence says otherwise.
              </p>
            </article>
            <div className="evidence-links" aria-label="Evidence and packaging links">
              {evidenceLinks.map(([label, href]) => (
                <a key={label} href={href}>{label}</a>
              ))}
            </div>
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
        <span>DotMatch</span>
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

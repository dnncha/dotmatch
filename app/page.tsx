const repoUrl = "https://github.com/dnncha/dotmatch";
const scientificClaimsUrl = `${repoUrl}/blob/main/docs/scientific-claims.md`;
const evidenceGalleryUrl = `${repoUrl}/blob/main/docs/evidence-gallery/README.md`;
const methodsUrl = `${repoUrl}/blob/main/docs/methods-and-citation.md`;
const packagingUrl = `${repoUrl}/blob/main/docs/packaging.md`;
const exposureUrl = `${repoUrl}/blob/main/docs/industry-exposure.md`;
const nextWinsUrl = `${repoUrl}/blob/main/docs/industry-next-wins.md`;
const workflowSubmissionsUrl = `${repoUrl}/blob/main/docs/workflow-submissions.md`;
const adoptersUrl = `${repoUrl}/blob/main/docs/adopters/README.md`;
const evaluationUrl = `${repoUrl}/blob/main/docs/bioinformatics-evaluation.md`;
const pilotUrl = `${repoUrl}/blob/main/docs/pilot-program.md`;
const reviewPacketUrl = `${repoUrl}/blob/main/docs/external-review-packet.md`;
const integrationRoadmapUrl = `${repoUrl}/blob/main/docs/workflow-integration-roadmap.md`;
const integrationTargetsUrl = `${repoUrl}/blob/main/docs/integration-targets.json`;
const reviewerReadinessUrl = `${repoUrl}/blob/main/docs/reviewer-readiness.json`;
const adoptionMetricsUrl = `${repoUrl}/blob/main/docs/adoption-metrics.md`;
const pypiUrl = "https://pypi.org/project/dotmatch/";

const basePath = process.env.NEXT_PUBLIC_BASE_PATH ?? "";
const assignmentWorkflowImage = `${basePath}/dotmatch-read-assignment.svg`;

const structuredData = {
  "@context": "https://schema.org",
  "@graph": [
    {
      "@type": "WebSite",
      "@id": "https://dnncha.github.io/dotmatch/#website",
      name: "AssayCode",
      alternateName: "DotMatch",
      url: "https://dnncha.github.io/dotmatch",
      description:
        "AssayCode compiles, validates, decodes, and diagnoses known-target sequencing assays using the DotMatch engine."
    },
    {
      "@type": "SoftwareApplication",
      "@id": "https://dnncha.github.io/dotmatch/#software",
      name: "AssayCode",
      alternateName: "DotMatch",
      applicationCategory: "Bioinformatics software",
      operatingSystem: "Linux, macOS",
      softwareHelp: "https://dotmatch.readthedocs.io/",
      codeRepository: repoUrl,
      license: `${repoUrl}/blob/main/LICENSE`,
      programmingLanguage: ["C", "Python", "R"],
      description:
        "AssayCode is an assay reliability platform powered by DotMatch for ambiguity-aware known-target DNA assignment."
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
  ["Claim boundaries", scientificClaimsUrl],
  ["Evidence gallery", evidenceGalleryUrl],
  ["Methods and citation", methodsUrl],
  ["Packaging notes", packagingUrl]
] as const;

const exposureActions = [
  {
    title: "Evaluate",
    body:
      "Install the released package, run the tutorial, and compare the explicit assignment outcomes against your current known-target workflow."
  },
  {
    title: "Integrate",
    body:
      "Use the workflow submission pack to make DotMatch visible in nf-core, MultiQC, Galaxy, Snakemake, and institutional pipeline reports."
  },
  {
    title: "Cite",
    body:
      "Copy methods language from the citation guidance so external reports describe assignment windows, ambiguity policy, and software version."
  },
  {
    title: "Pilot",
    body:
      "Record approved external pilots only after a public lab, workflow, or package integration can be linked and reviewed."
  },
  {
    title: "Share",
    body:
      "Use the outreach and integration kit for conference abstracts, repository announcements, short social copy, and direct maintainer outreach."
  }
] as const;

const nextWins = [
  "Decision tree",
  "Persona one-pagers",
  "Integration tracker",
  "Reviewer packet",
  "Conference abstracts",
  "Social pack",
  "Maintainer templates",
  "Pilot scorecard",
  "Adoption KPIs",
  "Release calendar"
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
          <a href="#industry-routes">Audiences</a>
          <a href="#evidence">Evidence</a>
          <a href="#exposure">Adoption</a>
          <a href="#next-wins">Next actions</a>
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

        <section id="exposure" className="section exposure-section" aria-labelledby="exposure-title">
          <div className="section-heading">
            <p className="section-kicker">Adoption flywheel</p>
            <h2 id="exposure-title">Five moves turn a useful tool into a visible one.</h2>
            <p>
              DotMatch is easier to recommend when every external mention points
              to runnable examples, scoped evidence, citation text, and a public
              record of accepted integrations.
            </p>
          </div>
          <div className="exposure-layout">
            <ol className="exposure-list">
              {exposureActions.map((action) => (
                <li key={action.title}>
                  <strong>{action.title}</strong>
                  <p>{action.body}</p>
                </li>
              ))}
            </ol>
            <aside className="exposure-links" aria-label="Adoption and exposure links">
          <a href={exposureUrl}>Outreach and integration kit</a>
              <a href={workflowSubmissionsUrl}>Workflow submission pack</a>
              <a href={evaluationUrl}>Bioinformatics evaluation packet</a>
              <a href={pilotUrl}>External pilot protocol</a>
              <a href={reviewPacketUrl}>External review packet</a>
              <a href={integrationRoadmapUrl}>Workflow integration roadmap</a>
              <a href={integrationTargetsUrl}>Integration target tracker</a>
              <a href={reviewerReadinessUrl}>Reviewer readiness record</a>
              <a href={adoptionMetricsUrl}>Adoption metrics contract</a>
              <a href={methodsUrl}>Methods and citation text</a>
              <a href={adoptersUrl}>Used-by record policy</a>
            </aside>
          </div>
        </section>

        <section id="next-wins" className="section next-wins-section" aria-labelledby="next-wins-title">
          <div className="section-heading">
            <p className="section-kicker">Next 10 exposure wins</p>
            <h2 id="next-wins-title">Turn interest into repeatable distribution.</h2>
            <p>
              The next adoption layer is a checked playbook: decision paths,
              persona-specific handoffs, maintainer-ready templates, pilot scoring,
              and release communication assets that can be reused without widening
              the scientific claim boundary.
            </p>
          </div>
          <div className="next-wins-layout">
            <ol className="next-wins-grid" aria-label="Next 10 distribution actions">
              {nextWins.map((win) => (
                <li key={win}>{win}</li>
              ))}
            </ol>
            <aside className="next-wins-note">
              <h3>Checked source of truth</h3>
              <p>
                The playbook and machine-readable tracker must stay aligned before
                release. Private outreach and unmerged PRs remain activity, not
                an accepted external-use record.
              </p>
              <a href={nextWinsUrl}>Open the next 10 playbook</a>
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

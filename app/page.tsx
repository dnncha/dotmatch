import { MobileNavigation } from "./mobile-navigation";
import packageMetadata from "../package.json";
import checkedFixtureEnvelope from "../agent-reference-crispr.json";

const repoUrl = "https://github.com/dnncha/dotmatch";
const docsUrl = "https://dotmatch.readthedocs.io/en/latest/";
const gettingStartedUrl = `${docsUrl}getting-started.html`;
const barcodeTroubleshootingUrl = `${gettingStartedUrl}#diagnose-a-barcode-run`;
const scopeUrl = `${docsUrl}trust-and-scope.html`;
const benchmarksUrl = `${docsUrl}benchmarks/README.html`;
const schemasUrl = `${docsUrl}schemas.html`;
const commandReferenceUrl = `${docsUrl}command-reference.html`;
const methodsUrl = `${docsUrl}methods-and-citation.html`;
const distributionUrl = `${docsUrl}release-process.html`;
const agentGuideUrl = `${docsUrl}agent-guide.html`;
const agentCrisprUrl = `${docsUrl}agent-crispr.html`;
const agentPerturbSeqUrl = `${docsUrl}agent-perturb-seq.html`;
const capabilityManifestUrl = "https://dnncha.github.io/dotmatch/agent-capabilities.json";
const agentToolsUrl = "https://dnncha.github.io/dotmatch/agent-tools.json";
const agentReferenceUrl = "https://dnncha.github.io/dotmatch/agent-reference-crispr.json";
const pypiUrl = "https://pypi.org/project/dotmatch/";
const containerUrl = `${repoUrl}/pkgs/container/dotmatch`;
const workflowExamplesUrl = `${repoUrl}/tree/main/examples/workflows`;
const doiUrl = "https://doi.org/10.5281/zenodo.20541628";
const binderUrl = "https://mybinder.org/v2/gh/dnncha/dotmatch/main?labpath=demo.ipynb";
const colabUrl = "https://colab.research.google.com/github/dnncha/dotmatch/blob/main/demo.ipynb";
const releaseVersion = packageMetadata.version;
const publishedVersion = "0.3.1";

const basePath = process.env.NEXT_PUBLIC_BASE_PATH ?? "";
const assignmentWorkflowImage = `${basePath}/dotmatch-read-assignment-v2.webp`;
const assignmentWorkflowMobileImage = `${basePath}/dotmatch-read-assignment-mobile-v2.webp`;

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
      softwareVersion: publishedVersion,
      softwareHelp: docsUrl,
      codeRepository: repoUrl,
      downloadUrl: pypiUrl,
      citation: doiUrl,
      license: `${repoUrl}/blob/main/LICENSE`,
      programmingLanguage: ["C", "Python", "R"],
      featureList: [
        "CRISPR guide counting",
        "inline barcode demultiplexing",
        "feature-barcode assignment",
        "CRISPR guide-capture assignment",
        "barcode panel design",
        "known-target FASTQ matching"
      ],
      subjectOf: capabilityManifestUrl,
      description:
        "DotMatch assigns fixed read windows to known short DNA targets and reports unique, ambiguous, none, and invalid outcomes for auditable sequencing workflows."
    }
  ]
};

const agentSteps = [
  ["01", "Discover", "Read the installed, versioned tool contract and scientific boundaries."],
  ["02", "Prepare", "Infer a fixed window and scaffold a reviewable local AssaySpec project."],
  ["03", "Preflight", "Audit target safety, resources, inference confidence, and the execution plan."],
  ["04", "Run", "Count locally, inspect reliability, and write a bounded numbered revision only when evidence supports one."],
  ["05", "Review + handoff", "Return normalized findings and a hashed bundle without copying raw FASTQ."]
] as const;

const agentTasks = [
  {
    id: "agent-crispr",
    title: "CRISPR guide counting",
    inputs: "Guide TSV/CSV, FASTQ directory, empty output directory",
    outputs: "MAGeCK-compatible counts, sample QC, reliability evidence, hashed handoff",
    limit: "Known-guide counting only; no downstream screen statistics.",
    start: `dotmatch agent invoke prepare_assay \\\n  --input crispr-request.json`,
    href: agentCrisprUrl
  },
  {
    id: "agent-perturb-seq",
    title: "Perturb-seq direct-guide capture",
    inputs: "Known guide-barcode table, guide-capture FASTQ directory, empty output directory",
    outputs: "Per-guide counts, per-read outcomes, reliability evidence, hashed handoff",
    limit: "No cell/UMI processing, guide-per-cell calls, expression, or perturbation effects.",
    start: `dotmatch agent invoke prepare_assay \\\n  --input perturb-seq-request.json`,
    href: agentPerturbSeqUrl
  }
] as const;

const outcomes = [
  ["unique", "One target fits. Count the read or write it to that target's output."],
  ["ambiguous", "Several targets fit. DotMatch makes no forced call."],
  ["none", "No target fits. Keep the read for unmatched-read review."],
  ["invalid", "The read window is missing. Record the failure in QC."]
] as const;

const failureModes = [
  {
    title: "Several targets fit",
    body:
      "Report the read as ambiguous instead of forcing it into one count."
  },
  {
    title: "The read window is missing",
    body:
      "Report short reads and invalid extraction windows instead of dropping them."
  },
  {
    title: "Correction could mix samples",
    body:
      "Check whether an error-correction setting could assign one sequence to several targets."
  },
  {
    title: "No expected target fits",
    body:
      "Keep recurring unmatched sequences available for assay, adapter, and off-target review."
  }
] as const;

const workflowSteps = [
  {
    step: "1",
    title: "Provide the sequences you expect",
    body:
      "Start with a table of guides, barcodes, feature tags, primers, panel targets, or whitelist sequences."
  },
  {
    step: "2",
    title: "Choose where to look in each read",
    body:
      "DotMatch extracts the same configured position from every read and compares it with the expected sequences."
  },
  {
    step: "3",
    title: "Review every outcome",
    body:
      "Inspect unique, ambiguous, unmatched, and invalid reads in counts, FASTQs, QC tables, and reports."
  }
] as const;

const audienceRoutes = [
  {
    title: "Core facilities",
    body:
      "Review ambiguous and unmatched sample barcodes, guide libraries, or panel targets before results leave the core.",
    link: "Start with barcode troubleshooting",
    href: barcodeTroubleshootingUrl
  },
  {
    title: "CRISPR screen teams",
    body:
      "Count known guide windows, write MAGeCK-compatible output, and keep assignment failures for methods review.",
    link: "Run the CRISPR tutorial",
    href: `${docsUrl}tutorials/crispr-count-first-run.html`
  },
  {
    title: "Workflow maintainers",
    body:
      "Carry stable TSV, JSON, FASTQ, and HTML output into nf-core, Galaxy, Snakemake, MultiQC, or local pipelines.",
    link: "See workflow examples",
    href: workflowExamplesUrl
  },
  {
    title: "Assay developers",
    body:
      "Design barcode panels, test correction radius safety, and export panel records before sequencing starts.",
    link: "Review panel design",
    href: `${docsUrl}barcode-panel-design.html`
  }
] as const;

const contexts = [
  "CRISPR guides",
  "inline barcodes",
  "feature tags",
  "primers / panels",
  "whitelists"
] as const;

const evidenceGroups = [
  {
    eyebrow: "Scope",
    title: "A narrow job, stated plainly",
    body:
      "DotMatch assigns a fixed window in each read to a known list of short DNA targets. It is not a genome aligner, basecaller, variant caller, UMI pipeline, or downstream CRISPR analysis package.",
    links: [
      ["Scope and limitations", scopeUrl],
      ["Methods and citation", methodsUrl]
    ]
  },
  {
    eyebrow: "Outputs",
    title: "Files you can inspect",
    body:
      "Review count tables, split FASTQs, QC tables, JSON summaries, and HTML reports. Stable schemas and workflow examples make the results easier to carry into a pipeline.",
    links: [
      ["Output schemas", schemasUrl],
      ["Workflow examples", workflowExamplesUrl]
    ]
  },
  {
    eyebrow: "Release",
    title: `Published ${publishedVersion}; candidate ${releaseVersion}`,
    body:
      "Version 0.3.1 is the current verified release on PyPI, GHCR, GitHub Releases, and Zenodo. Version 0.4.0 is a source candidate until it is separately authorized and published; Bioconda and BioContainers remain downstream gates.",
    links: [
      ["Release status", distributionUrl],
      ["Benchmarks", benchmarksUrl]
    ]
  }
] as const;

const trustFacts = [
  ["Candidate", `Version ${releaseVersion} source candidate`],
  ["License", "Apache-2.0 open source"],
  ["Data", "Runs locally on your machine"],
  ["Citation", "Archived release with DOI"]
] as const;

const navigationLinks = [
  ["Run with an agent", "#agent-workflow"],
  ["How it works", "#workflow"],
  ["Reliability", "#failure-modes"],
  ["Uses", "#use-cases"],
  ["Methods", "#evidence"],
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
            <p className="positioning">Known-target sequencing read assignment</p>
            <h1 id="hero-title">Match reads without hiding uncertainty.</h1>
            <p className="hero-lede">
              Compare the same fixed window in every read with the short DNA sequences
              you expect. DotMatch labels each read unique, ambiguous, unmatched, or invalid.
            </p>
            <p className="hero-text">
              Use it for CRISPR guides, inline barcodes, feature tags, primers,
              panel targets, and whitelist-style assays. It runs locally. It is not
              a genome aligner or basecaller.
            </p>
            <div className="hero-actions">
              <a className="button primary" href="#agent-workflow">Run with a local agent</a>
              <a className="button secondary" href="#install">Install the CLI</a>
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
              <picture>
                <source
                  media="(max-width: 520px)"
                  srcSet={assignmentWorkflowMobileImage}
                  width="864"
                  height="1821"
                />
                <img
                  src={assignmentWorkflowImage}
                  width="1825"
                  height="862"
                  alt="DotMatch takes FASTQ reads and a known target list, compares the same fixed-position window in every read, and reports each read as unique, ambiguous, unmatched, or invalid"
                  decoding="async"
                  fetchPriority="high"
                />
              </picture>
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

        <section id="agent-workflow" className="section agent-section" aria-labelledby="agent-title">
          <div className="section-heading">
            <p className="section-kicker">A local research agent, with evidence</p>
            <h2 id="agent-title">From assay files to a reviewable handoff.</h2>
            <p>
              A Codex or other local agent can use six structured DotMatch tools to
              prepare, run, correct, and package a known-guide workflow. It receives
              paths and scientific parameters—not shell commands—and ordinary runs
              make no network requests.
            </p>
          </div>
          <ol className="agent-steps" aria-label="Local agent workflow">
            {agentSteps.map(([step, title, body]) => (
              <li key={step}>
                <span>{step}</span>
                <strong>{title}</strong>
                <p>{body}</p>
              </li>
            ))}
          </ol>
          <div className="agent-start-grid">
            <article className="agent-command">
              <h3>Start from the installed contract</h3>
              <pre><code>{`dotmatch agent tools --json
dotmatch agent export-skill --target ./dotmatch-agent`}</code></pre>
              <div className="agent-links">
                <a href={agentGuideUrl}>Agent guide</a>
                <a href={agentToolsUrl}>Tool contract JSON</a>
                <a href={capabilityManifestUrl}>Capability manifest 1.1</a>
              </div>
            </article>
            <article className="agent-evidence">
              <h3>Checked executable fixture, honest verdict</h3>
              <p>
                This checked CRISPR fixture contains unique, ambiguous, unmatched,
                and invalid reads. Its expected counts and sample-QC artifacts are
                hashed; the structured verdict remains <code>failed</code>
                because the fixture intentionally exercises unsafe and low-assignment states.
              </p>
              <pre><code>{JSON.stringify(checkedFixtureEnvelope, null, 2)}</code></pre>
              <a href={agentReferenceUrl}>Open the checked fixture envelope</a>
            </article>
          </div>
          <div className="agent-task-grid" aria-label="Supported agent task pages">
            {agentTasks.map((task) => (
              <article id={task.id} key={task.id}>
                <h3>{task.title}</h3>
                <dl>
                  <div><dt>Inputs</dt><dd>{task.inputs}</dd></div>
                  <div><dt>Outputs</dt><dd>{task.outputs}</dd></div>
                  <div><dt>Limit</dt><dd>{task.limit}</dd></div>
                </dl>
                <pre><code>{task.start}</code></pre>
                <a href={task.href}>Open the exact task page</a>
              </article>
            ))}
          </div>
        </section>

        <section id="failure-modes" className="section failure-section" aria-labelledby="failure-title">
          <div className="section-heading">
            <p className="section-kicker">No silent failures</p>
            <h2 id="failure-title">Keep the hard calls visible.</h2>
            <p>
              A target list does not make every read safe to count. DotMatch separates
              uncertain and invalid reads so you can review them instead of losing them.
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
            <h2 id="workflow-title">From target list to auditable result in three steps.</h2>
            <p>
              Provide the sequences you expect, choose where to look in each read,
              then inspect every assignment.
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

        <section id="use-cases" className="section audience-section" aria-labelledby="audience-title">
          <div className="section-heading">
            <p className="section-kicker">Choose your path</p>
            <h2 id="audience-title">Start with the work in front of you.</h2>
            <p>
              DotMatch fits assays where the short target sequences are known in
              advance and every assignment needs to remain inspectable.
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
            <p className="section-kicker">Methods, outputs, and release</p>
            <h2 id="evidence-title">Know exactly what DotMatch does—and does not do.</h2>
            <p>
              The useful details are close at hand: scientific scope, inspectable
              outputs, benchmark methods, and the status of each release channel.
            </p>
          </div>
          <div className="evidence-grid">
            {evidenceGroups.map((group) => (
              <article key={group.eyebrow} className="evidence-card">
                <p className="evidence-eyebrow">{group.eyebrow}</p>
                <h3>{group.title}</h3>
                <p>{group.body}</p>
                <div className="evidence-card-links">
                  {group.links.map(([label, href]) => (
                    <a key={label} href={href}>{label}</a>
                  ))}
                </div>
              </article>
            ))}
          </div>
          <p className="automation-note">
            Building an automated workflow? See the <a href={agentGuideUrl}>agent task guide</a>,
            {" "}the <a href={agentToolsUrl}>versioned tool contract</a>, the
            {" "}<a href={capabilityManifestUrl}>machine-readable capabilities</a>, or the
            {" "}<a href={commandReferenceUrl}>command reference</a>.
          </p>
          <p className="automation-note">
            Snakemake wrapper PR #5825 is accepted on its default branch but absent from
            v9.16.0. The nf-core, Galaxy, and MultiQC contributions remain under review.
          </p>
        </section>

        <section id="install" className="section install-section" aria-labelledby="install-title">
          <div className="install-copy">
            <p className="section-kicker">Try it locally</p>
            <h2 id="install-title">Start with a two-sequence check.</h2>
            <p>
              Install the command-line package, then compare two short sequences
              before you point DotMatch at FASTQ data. Your files stay on your machine.
            </p>
            <p>
              Prefer to try it in a browser? Launch the small synthetic DotMatch
              notebook in <a href={binderUrl}>Binder</a> or <a href={colabUrl}>Google
              Colab</a>. No local setup is required.
            </p>
          </div>
          <div className="terminal" aria-label="Install DotMatch">
            <div className="terminal-bar" aria-hidden="true">
              <span />
              <span />
              <span />
            </div>
            <pre><code>{`python3 -m pip install dotmatch
dotmatch dist ACGT AGGT`}</code></pre>
            <div className="terminal-links">
              <a href={pypiUrl}>PyPI package</a>
              <a href={gettingStartedUrl}>Getting started</a>
              <a href={containerUrl}>Pinned container</a>
              <a href={distributionUrl}>Release status</a>
              <a href={methodsUrl}>Citation guidance</a>
              <a href={binderUrl}>Launch Binder demo</a>
              <a href={colabUrl}>Launch Colab demo</a>
            </div>
          </div>
        </section>
      </main>

      <footer className="site-footer">
        <span>DotMatch · Apache-2.0</span>
        <nav aria-label="Footer navigation">
          <a href={docsUrl}>Documentation</a>
          <a href={repoUrl}>GitHub</a>
          <a href={scopeUrl}>Scope and limitations</a>
          <a href={methodsUrl}>Methods</a>
          <a href={doiUrl}>DOI</a>
        </nav>
      </footer>
    </>
  );
}

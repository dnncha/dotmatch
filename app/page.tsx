import packageMetadata from "../package.json";
import demo from "../public/assignment-demo.json";
import { ResearchHeader, ResearchFooter } from "./research-shell";
import {
  canonicalUrl,
  conceptDoi,
  docsUrl,
  publishedVersion,
  repoUrl,
  sitePath,
  siteAsset,
} from "./site-metadata";
import { AssignmentDemo } from "./assignment-demo";
import { InstallCommand } from "./install-command";
import { GuideCounterBenchmark } from "./guide-counter-benchmark";
import styles from "./home.module.css";

const structuredData = {
  "@context": "https://schema.org",
  "@graph": [
    {
      "@type": "WebSite",
      "@id": `${canonicalUrl()}#website`,
      name: "DotMatch",
      url: canonicalUrl(),
      description:
        "CRISPR guide counting and barcode assignment with explicit read outcomes.",
    },
    {
      "@type": "SoftwareApplication",
      "@id": `${canonicalUrl()}#software`,
      name: "DotMatch",
      applicationCategory: "Bioinformatics software",
      operatingSystem: "Linux, macOS",
      softwareVersion: publishedVersion,
      subjectOf: `${canonicalUrl()}agent-capabilities.json`,
      softwareHelp: docsUrl,
      codeRepository: repoUrl,
      downloadUrl: "https://pypi.org/project/dotmatch/",
      citation: conceptDoi,
      license: `${repoUrl}/blob/main/LICENSE`,
      programmingLanguage: ["C", "Python", "R"],
      featureList: [
        "CRISPR guide counting",
        "MAGeCK-compatible count tables",
        "inline barcode demultiplexing",
        "target-library auditing",
        "explicit ambiguous-read outcomes",
      ],
      description:
        "Assign fixed read windows from FASTQ to known CRISPR guides, barcodes, and short DNA targets while keeping unique, ambiguous, unmatched, and invalid outcomes visible.",
    },
  ],
};

const policyRows = [
  ["exact", "Exact"],
  ["radius_k1", "One target within k=1"],
  ["best_k1", "Nearest target, k=1"],
] as const;

export default function Home() {
  return (
    <>
      <script
        type="application/ld+json"
        dangerouslySetInnerHTML={{ __html: JSON.stringify(structuredData) }}
      />
      <ResearchHeader />
      <main id="main-content" className={styles.home}>
        <section id="top" className={styles.hero} aria-labelledby="hero-title">
          <div>
            <h1 id="hero-title">
              CRISPR guide counting and barcode assignment you can audit
            </h1>
            <p className={styles.lede}>
              DotMatch matches a fixed window from every FASTQ read against your
              known target library locally on Linux and macOS. Export
              MAGeCK-ready count tables or split FASTQs while keeping ambiguous,
              unmatched, and invalid reads visible.
            </p>
            <div className={styles.actions}>
              <a
                className={styles.primary}
                href={sitePath("crispr-guide-counting")}
              >
                Count CRISPR guides
              </a>
              <a
                className={styles.secondary}
                href={sitePath("tools/library-safety")}
              >
                Check a library
              </a>
            </div>
            <div className={styles.secondaryLinks}>
              <a
                className={styles.textLink}
                href={siteAsset("examples/assignment-review/report.html")}
              >
                Open example report
              </a>
              <a
                className={styles.textLink}
                href={`${docsUrl}tutorials/crispr-count-first-run.html`}
              >
                First-run tutorial
              </a>
            </div>
            <p className={styles.facts}>
              <span>Apache-2.0</span>
              <span>Linux and macOS</span>
              <span>Runs locally (no uploads)</span>
            </p>
          </div>
          <AssignmentDemo />
        </section>

        <section
          id="install"
          className={styles.installStrip}
          aria-labelledby="install-title"
        >
          <div>
            <h2 id="install-title">Install DotMatch</h2>
            <p>
              Published CLI {publishedVersion} ·{" "}
              <a href={`${docsUrl}getting-started.html`}>
                Conda, containers, and source builds
              </a>
            </p>
          </div>
          <InstallCommand />
        </section>

        <GuideCounterBenchmark />

        <section
          id="failure-modes"
          className={styles.section}
          aria-labelledby="sensitivity-title"
        >
          <div className={styles.proof}>
            <div>
              <h2 id="sensitivity-title">
                Why mapped-read totals hide count shifts
              </h2>
              <p>
                In this nine-read synthetic fixture, exact matching and
                radius-one matching both produce three unique assignments, yet
                they disagree on individual guide counts.
              </p>
              <div className={styles.proofLinks}>
                <a
                  className={styles.textLink}
                  href={siteAsset("examples/assignment-review/report.html")}
                >
                  Open example report
                </a>
                <a
                  className={styles.textLink}
                  href={sitePath("assignment-sensitivity")}
                >
                  Explore sensitivity
                </a>
              </div>
            </div>
            <div>
              <div
                className={styles.proofTable}
                tabIndex={0}
                role="region"
                aria-label="Checked synthetic policy comparison"
              >
                <table>
                  <caption>
                    Synthetic fixture ({demo.read_count} reads,{" "}
                    {demo.target_count} target IDs)
                  </caption>
                  <thead>
                    <tr>
                      <th scope="col">Matching rule</th>
                      <th scope="col">Unique</th>
                      <th scope="col">Ambiguous</th>
                      <th scope="col">Unmatched</th>
                    </tr>
                  </thead>
                  <tbody>
                    {policyRows.map(([key, label]) => (
                      <tr key={key}>
                        <th scope="row">{label}</th>
                        <td>{demo.outcomes[key].unique}</td>
                        <td>{demo.outcomes[key].ambiguous}</td>
                        <td>{demo.outcomes[key].none}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              <p className={styles.proofNote}>
                One additional read has an invalid window under every policy, and{" "}
                {demo.changed_reads} reads change outcome across policies.
                Calculated with the native matcher to demonstrate assignment
                behavior rather than biological accuracy.
              </p>
            </div>
          </div>
        </section>

        <section
          id="workflow"
          className={styles.section}
          aria-labelledby="workflow-title"
        >
          <h2 id="workflow-title">How a counting run works</h2>
          <ol className={styles.steps}>
            <li>
              <span className={styles.stepMarker} aria-hidden="true">
                1
              </span>
              <h3>Confirm library and read window</h3>
              <p>
                Provide your guide or barcode table alongside FASTQ files, and
                review the inferred window and orientation before counting.
              </p>
              <a href={`${docsUrl}tutorials/crispr-count-first-run.html`}>
                Follow the first-run tutorial
              </a>
            </li>
            <li>
              <span className={styles.stepMarker} aria-hidden="true">
                2
              </span>
              <h3>Assign reads under explicit rules</h3>
              <p>
                Run exact, Hamming, or Levenshtein matching while keeping
                ambiguous and unmatched reads out of target counts.
              </p>
              <a href={sitePath("assignment-sensitivity")}>
                Compare matching policies
              </a>
            </li>
            <li>
              <span className={styles.stepMarker} aria-hidden="true">
                3
              </span>
              <h3>Export counts and a portable QC bundle</h3>
              <p>
                Pass MAGeCK-ready counts to downstream statistics, and archive
                config, QC, and SHA-256 provenance without copying raw FASTQs.
              </p>
              <a href={`${docsUrl}lab-evaluation.html`}>
                Inspect the handoff workflow
              </a>
            </li>
          </ol>
        </section>

        <section
          id="use-cases"
          className={styles.section}
          aria-labelledby="use-title"
        >
          <h2 id="use-title">Where DotMatch fits</h2>
          <div className={styles.routes}>
            <article>
              <h3>Pooled CRISPR screens</h3>
              <p>
                Count known guides, review assignment QC, and export a
                MAGeCK-compatible guide-by-sample matrix.
              </p>
              <a href={sitePath("crispr-guide-counting")}>CRISPR counting</a>
            </article>
            <article>
              <h3>Inline barcodes and sequencing cores</h3>
              <p>
                Split pooled FASTQs by inline barcode and run{" "}
                <code>barcode autopsy</code> to diagnose shifted offsets or
                unsafe correction.
              </p>
              <a href={`${docsUrl}getting-started.html#diagnose-a-barcode-run`}>
                Barcode diagnostics
              </a>
            </article>
            <article>
              <h3>Feature and guide capture</h3>
              <p>
                Assign known feature tags and build sparse cell-by-feature
                matrices from extracted observations with cell IDs.
              </p>
              <a href={`${docsUrl}tutorials/scverse-perturb-seq.html`}>
                Feature-matrix workflow
              </a>
            </article>
          </div>
        </section>

        <section
          id="evidence"
          className={styles.section}
          aria-labelledby="evidence-title"
        >
          <h2 id="evidence-title">Benchmarks and reproducible examples</h2>
          <p>
            Reproduce public comparisons, inspect count differences, and check
            DotMatch alongside your existing pipeline.
          </p>
          <div className={styles.evidenceRows}>
            <article className={styles.evidenceRow}>
              <h3>Yusa and Brunello CRISPR data</h3>
              <p>
                Paired comparisons against MAGeCK, guide-counter, and reference
                matchers reporting runtime, memory, parameters, and count
                differences together.
              </p>
              <a href={`${docsUrl}benchmarks/crispr_comparison/README.html`}>
                Read the results
              </a>
            </article>
            <article className={styles.evidenceRow}>
              <h3>GSE146194 direct-guide capture</h3>
              <p>
                Public benchmark separating discovery and evaluation reads to
                check per-read guide assignment (separate from full single-cell
                analysis).
              </p>
              <a href={`${repoUrl}/tree/main/examples/perturb_seq_gse146194`}>
                Reproduce the case study
              </a>
            </article>
            <article className={styles.evidenceRow}>
              <h3>Nine-read sensitivity fixture</h3>
              <p>
                Checked fixture exercising close targets, duplicate IDs, a
                literal N, an unmatched read, and a short read against the
                native matcher.
              </p>
              <a href={`${repoUrl}/tree/main/examples/assignment_sensitivity`}>
                Inspect the fixture
              </a>
            </article>
          </div>
        </section>

        <section
          id="agent-workflow"
          className={styles.automation}
          aria-label="Automation and scientific scope"
        >
          <details>
            <summary>Pipelines, MCP server, and local AI agents</summary>
            <p>
              Connect DotMatch as a local MCP server with{" "}
              <code>dotmatch mcp</code>, or inspect the structured CLI tools and
              output schemas before wiring a pipeline.
            </p>
            <pre>
              <code>{`dotmatch agent tools --json\ndotmatch agent export-skill --target ./dotmatch-agent`}</code>
            </pre>
            <p>
              <a href={`${docsUrl}agent-guide.html`}>Agent guide</a> ·{" "}
              <a href={`${repoUrl}/tree/main/examples/workflows`}>
                Workflow examples
              </a>{" "}
              ·{" "}
              <a href={`${docsUrl}command-reference.html`}>Command reference</a>
            </p>
          </details>
          <p className={styles.boundary}>
            DotMatch handles fixed-window known-target sequence assignment and
            powers{" "}
            <a href={`${docsUrl}crisprworks.html`}>CRISPRWorks Count</a>. Genome
            alignment, cell/UMI processing, and downstream screen statistics
            remain separate steps.{" "}
            <a href={`${docsUrl}trust-and-scope.html`}>Methods and scope</a> ·{" "}
            <a href={`${docsUrl}methods-and-citation.html`}>Cite DotMatch</a>
            <span className={styles.boundaryVersion}>
              Published package {publishedVersion}; website source version{" "}
              {packageMetadata.version}.
            </span>
          </p>
        </section>
      </main>
      <ResearchFooter />
    </>
  );
}

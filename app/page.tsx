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
      name: "CRISPRWorks Count, powered by DotMatch",
      url: canonicalUrl(),
      description:
        "Count CRISPR guides from FASTQ files, export MAGeCK-compatible tables, and review how reads were assigned.",
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
        "Count known CRISPR guides from FASTQ and inspect unique, ambiguous, unmatched and invalid read outcomes before downstream analysis.",
    },
  ],
};
const policyRows = [
  ["exact", "Exact match"],
  ["radius_k1", "One target within 1 substitution"],
  ["best_k1", "Unique nearest target within 1 substitution"],
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
            <p className={styles.eyebrow}>
              CRISPRWorks Count, powered by DotMatch
            </p>
            <h1 id="hero-title">
              Count CRISPR guides from FASTQ files
            </h1>
            <p className={styles.lede}>
              Turn your guide library and sequencing reads into MAGeCK-compatible
              count tables. Review which reads count toward a guide and which
              remain ambiguous, unmatched, or invalid.
            </p>
            <div className={styles.actions}>
              <a
                className={styles.primary}
                href={sitePath("crispr-guide-counting")}
              >
                Start counting guides
              </a>
              <a
                className={styles.textLink}
                href={sitePath("tools/library-safety")}
              >
                Check library collisions
              </a>
              <a
                className={styles.textLink}
                href={`${docsUrl}tutorials/crispr-count-first-run.html`}
              >
                Run the synthetic example
              </a>
              <a
                className={styles.textLink}
                href={siteAsset("examples/assignment-review/report.html")}
              >
                Open example report
              </a>
            </div>
            <p className={styles.facts}>
              <span>Apache-2.0</span>
              <span>Linux & macOS</span>
              <span>Local processing</span>
            </p>
            <p>
              Install and run as <code>dotmatch</code>. Part of the{" "}
              <a href={`${docsUrl}crisprworks.html`}>CRISPRWorks family</a>.
            </p>
          </div>
          <AssignmentDemo />
        </section>
        <GuideCounterBenchmark />
        <section
          id="install"
          className={styles.installStrip}
          aria-labelledby="install-title"
        >
          <div>
            <h2 id="install-title">Install the command-line tool</h2>
            <p>
              DotMatch {publishedVersion} ·{" "}
              <a href={`${docsUrl}getting-started.html`}>
                Conda, containers, and installation help
              </a>
            </p>
          </div>
          <InstallCommand />
        </section>
        <section
          id="workflow"
          className={styles.section}
          aria-labelledby="workflow-title"
        >
          <h2 id="workflow-title">From reads to count tables</h2>
          <ol className={styles.steps}>
            <li>
              <span className={styles.stepNumber}>01 / PREPARE</span>
              <h3>Prepare inputs</h3>
              <p>
                Provide a guide library and the FASTQ files containing your guide
                reads. Confirm the proposed guide window, orientation, and
                sample mapping before counting.
              </p>
              <a href={`${docsUrl}tutorials/crispr-count-first-run.html`}>
                Follow the first-run tutorial
              </a>
            </li>
            <li>
              <span className={styles.stepNumber}>02 / COUNT</span>
              <h3>Count matching reads</h3>
              <p>
                Choose exact matching or allow substitutions or insertions and
                deletions. Only reads assigned to one target add to its count.
              </p>
              <a href={sitePath("assignment-sensitivity")}>
                Compare matching policies
              </a>
            </li>
            <li>
              <span className={styles.stepNumber}>03 / REVIEW</span>
              <h3>Review counts and quality checks</h3>
              <p>
                Inspect read outcomes and guide counts. Export the raw count
                matrix for MAGeCK and keep the configuration, quality checks,
                and methods together in a local review bundle.
              </p>
              <a href={`${docsUrl}lab-evaluation.html`}>
                Review and share a completed run
              </a>
            </li>
          </ol>
        </section>
        <section
          id="failure-modes"
          className={styles.section}
          aria-labelledby="sensitivity-title"
        >
          <div className={styles.proof}>
            <div>
              <p className={styles.eyebrow}>
                Matching rules and guide counts
              </p>
              <h2 id="sensitivity-title">
                Check what changes when you allow a mismatch
              </h2>
              <p>
                In this synthetic nine-read example, exact matching and a rule
                allowing one substitution each assign three reads. The counts
                differ for three guides. The assigned total alone hides that
                change.
              </p>
              <a
                className={styles.textLink}
                href={siteAsset("examples/assignment-review/report.html")}
              >
                Open example report
              </a>
            </div>
            <div>
              <div
                className={styles.proofTable}
                tabIndex={0}
                role="region"
                aria-label="Matching rules compared on synthetic reads"
              >
                <table>
                  <caption>
                    Synthetic example · {demo.read_count} reads ·{" "}
                    {demo.target_count} target IDs
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
                One additional read is too short for the required window under
                every rule. {demo.changed_reads} reads change outcome between
                rules. These results come from the native matcher and demonstrate
                software behavior. They do not measure biological accuracy.
              </p>
            </div>
          </div>
        </section>
        <section
          id="use-cases"
          className={styles.section}
          aria-labelledby="use-title"
        >
          <h2 id="use-title">Choose your task</h2>
          <div className={styles.routes}>
            <article>
              <h3>Pooled CRISPR screens</h3>
              <p>
                Count known guides, review read assignments, and export a
                guide-by-sample matrix.
              </p>
              <a href={sitePath("crispr-guide-counting")}>Count CRISPR guides</a>
            </article>
            <article>
              <h3>Inline barcode demultiplexing</h3>
              <p>
                Demultiplex inline barcodes and investigate collisions or
                unexpected unmatched reads.
              </p>
              <a href={`${docsUrl}getting-started.html#diagnose-a-barcode-run`}>
                Investigate a barcode run
              </a>
            </article>
            <article>
              <h3>Feature and guide capture</h3>
              <p>
                Build cell-by-feature matrices from extracted feature windows
                and existing cell identifiers. Cell calling and UMI
                deduplication are separate upstream steps.
              </p>
              <a href={`${docsUrl}tutorials/scverse-perturb-seq.html`}>
                Build a feature matrix
              </a>
            </article>
          </div>
        </section>
        <section
          id="evidence"
          className={styles.section}
          aria-labelledby="evidence-title"
        >
          <h2 id="evidence-title">Benchmarks and examples</h2>
          <p>
            Reproduce a public example, inspect count differences, and compare
            DotMatch with your current workflow.
          </p>
          <div className={styles.evidenceRows}>
            <article className={styles.evidenceRow}>
              <h3>Yusa and Brunello CRISPR data</h3>
              <p>
                Recorded comparisons with MAGeCK, guide-counter and reference
                matchers. Runtime, memory, settings and count differences are
                reported together.
              </p>
              <a href={`${docsUrl}benchmarks/crispr_comparison/README.html`}>
                Read the CRISPR comparison results
              </a>
            </article>
            <article className={styles.evidenceRow}>
              <h3>Public direct-guide capture</h3>
              <p>
                A GSE146194 example with separate discovery and evaluation
                reads. Evidence for per-read guide assignment, not completed
                single-cell analysis.
              </p>
              <a href={`${repoUrl}/tree/main/examples/perturb_seq_gse146194`}>
                Reproduce the guide-capture case study
              </a>
            </article>
            <article className={styles.evidenceRow}>
              <h3>The example on this page</h3>
              <p>
                Nine synthetic reads exercise close targets, duplicate
                sequences, a literal N, an unmatched read and a short read.
                Expected assignments are checked against the native matcher.
              </p>
              <a href={`${repoUrl}/tree/main/examples/assignment_sensitivity`}>
                Inspect the synthetic reads and expected counts
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
            <summary>Pipelines and local agents</summary>
            <p>
              Prepare an assay, check its settings, run it, and review the
              results with local agent tools. List the tools supported by your
              installed version before choosing a workflow.
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
            DotMatch assigns reads to known short-DNA targets. Use separate tools
            for genome alignment, cell calling, UMI deduplication, and gene-level
            screen analysis.{" "}
            <a href={`${docsUrl}trust-and-scope.html`}>Methods & scope</a> ·{" "}
            <a href={`${docsUrl}methods-and-citation.html`}>Cite DotMatch</a>
            <br />
            Published package {publishedVersion}; website source version{" "}
            {packageMetadata.version}.
          </p>
        </section>
      </main>
      <ResearchFooter />
    </>
  );
}

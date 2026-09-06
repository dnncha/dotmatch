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
} from "./site-metadata";
import { AssignmentDemo } from "./assignment-demo";
import { InstallCommand } from "./install-command";
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
        "Count known CRISPR guides from FASTQ and inspect unique, ambiguous, unmatched and invalid read outcomes before downstream analysis.",
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
            <p className={styles.eyebrow}>
              CRISPR guide counting & barcode assignment
            </p>
            <h1 id="hero-title">
              Count your guides.<span>Account for every read.</span>
            </h1>
            <p className={styles.lede}>
              Go from FASTQ and a guide library to MAGeCK-compatible counts. See
              which reads matched, which were ambiguous, and which need another
              look.
            </p>
            <div className={styles.actions}>
              <a
                className={styles.primary}
                href={sitePath("crispr-guide-counting")}
              >
                Start counting guides <span aria-hidden="true">&nbsp;→</span>
              </a>
              <a
                className={styles.textLink}
                href={sitePath("tools/library-safety")}
              >
                Check a library first
              </a>
            </div>
            <p className={styles.facts}>
              <span>Apache-2.0</span>
              <span>Linux & macOS</span>
              <span>Your data stays local</span>
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
            <h2 id="install-title">Install. Bring your own data.</h2>
            <p>
              Published CLI {publishedVersion} ·{" "}
              <a href={`${docsUrl}getting-started.html`}>
                Conda, containers & installation help
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
          <p className={styles.eyebrow}>From reads to a count matrix</p>
          <h2 id="workflow-title">
            A small step in your workflow.
            <br />A clear account of your data.
          </h2>
          <ol className={styles.steps}>
            <li>
              <span className={styles.stepNumber}>01 / PREPARE</span>
              <h3>Start with the guides you expect.</h3>
              <p>
                Give DotMatch your library and FASTQs. Review the proposed read
                window before the analysis starts.
              </p>
              <a href={`${docsUrl}tutorials/crispr-count-first-run.html`}>
                Follow the first-run tutorial
              </a>
            </li>
            <li>
              <span className={styles.stepNumber}>02 / COUNT</span>
              <h3>Make the matching rule explicit.</h3>
              <p>
                Use exact, substitution-tolerant or indel-aware matching. Keep
                ambiguous and unmatched reads separate from unique counts.
              </p>
              <a href={sitePath("assignment-sensitivity")}>
                See why the rule matters
              </a>
            </li>
            <li>
              <span className={styles.stepNumber}>03 / ANALYSE</span>
              <h3>Keep the analysis you already trust.</h3>
              <p>
                Take the raw count matrix into MAGeCK. Carry configuration, QC
                and methods with the result in a local review bundle.
              </p>
              <a href={`${docsUrl}lab-evaluation.html`}>
                Inspect the handoff workflow
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
                More assigned reads ≠ better assignments
              </p>
              <h2 id="sensitivity-title">
                Same reads.
                <br />
                Different rules.
                <br />
                Different counts.
              </h2>
              <p>
                Here, exact matching and radius-one matching both count three
                reads—but not the same reads or guides. A single mapping
                percentage would hide that difference.
              </p>
              <a
                className={styles.textLink}
                href={sitePath("assignment-sensitivity")}
              >
                Explore the worked example →
              </a>
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
                One additional read has an invalid window under every policy.{" "}
                {demo.changed_reads} reads change outcome between policies.
                Calculated with the native matcher; not a biological accuracy
                benchmark.
              </p>
            </div>
          </div>
        </section>
        <section
          id="use-cases"
          className={styles.section}
          aria-labelledby="use-title"
        >
          <h2 id="use-title">Built around known sequences.</h2>
          <div className={styles.routes}>
            <article>
              <h3>Pooled CRISPR screens</h3>
              <p>
                Count known guides, inspect assignment QC and export a
                guide-by-sample matrix.
              </p>
              <a href={sitePath("crispr-guide-counting")}>CRISPR counting</a>
            </article>
            <article>
              <h3>Barcodes & sequencing cores</h3>
              <p>
                Demultiplex inline barcodes and investigate collisions or
                unexpected unmatched reads.
              </p>
              <a href={`${docsUrl}getting-started.html#diagnose-a-barcode-run`}>
                Barcode diagnostics
              </a>
            </article>
            <article>
              <h3>Feature & guide capture</h3>
              <p>
                Assign known features and build matrices from observations with
                explicit cell identifiers.
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
          <p className={styles.eyebrow}>Examine it before you adopt it</p>
          <h2 id="evidence-title">The inputs, methods and results are open.</h2>
          <p>
            Reproduce a public example, inspect count differences and compare
            DotMatch with your current workflow.
          </p>
          <div className={styles.evidenceRows}>
            <article className={styles.evidenceRow}>
              <h3>Yusa & Brunello CRISPR data</h3>
              <p>
                Recorded comparisons with MAGeCK, guide-counter and reference
                matchers. Runtime, memory, settings and count differences are
                reported together.
              </p>
              <a href={`${docsUrl}benchmarks/crispr_comparison/README.html`}>
                Read the results →
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
                Reproduce the case study →
              </a>
            </article>
            <article className={styles.evidenceRow}>
              <h3>The example on this page</h3>
              <p>
                Nine synthetic reads exercise close targets, duplicate
                sequences, a literal N, an unmatched read and a short read.
                Every displayed result is checked from source.
              </p>
              <a href={`${repoUrl}/tree/main/examples/assignment_sensitivity`}>
                Inspect the fixture →
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
            <summary>Working in a pipeline or with a local agent?</summary>
            <p>
              Use stable files and structured tools to prepare, preflight, run
              and review an assay. Inspect the installed contract before
              choosing a workflow.
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
            DotMatch handles known-target sequence assignment. Genome alignment,
            cell/UMI processing and downstream screen statistics remain separate
            steps.{" "}
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

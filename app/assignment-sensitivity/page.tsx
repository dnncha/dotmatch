import demo from "../../public/assignment-demo.json";
import { ResearchHeader, ResearchFooter } from "../research-shell";
import {
  docsUrl,
  repoUrl,
  sitePath,
  siteAsset,
  pageMetadata,
  publishedVersion,
} from "../site-metadata";
import styles from "../research.module.css";
export const metadata = pageMetadata(
  "Assignment sensitivity | DotMatch",
  "Compare exact, radius-one Hamming and best-distance Hamming assignments. View guide counts, read transitions and source artifacts.",
  "assignment-sensitivity",
);
const modes = [
  ["exact", "Exact"],
  ["radius_k1", "Radius one"],
  ["best_k1", "Best distance, k=1"],
] as const;
export default function SensitivityPage() {
  return (
    <>
      <ResearchHeader />
      <main id="main-content" className={styles.shell}>
        <header className={styles.pageHeader}>
          <p className={styles.breadcrumb}>
            <a href={sitePath()}>DotMatch</a> / Assignment sensitivity
          </p>
          <h1 className={styles.title}>
            Assignment sensitivity
          </h1>
          <p className={styles.lede}>
            Compare exact, radius-one Hamming and best-distance Hamming
            assignments on a nine-read synthetic dataset.
          </p>
          <div className={styles.actions}>
            <a className={styles.primary} href={siteAsset("examples/assignment-review/report.html")}>
              Open example report
            </a>
            <a className={styles.secondary} href={siteAsset("examples/assignment-review/dotmatch-review-example.zip")} download>
              Download example files
            </a>
          </div>
          <p className={styles.note}>
            The report uses results computed by the native matcher. The site
            viewer reviews the bundled synthetic example. The installed
            sensitivity command still writes its static report.
          </p>
        </header>
        <section className={styles.section} aria-labelledby="try-review">
          <h2 id="try-review">Using the example report</h2>
          <div className={styles.grid}>
            <article className={styles.card}>
              <h3>1. Compare counts</h3>
              <p>Exact and Radius k=1 both assign three reads. Their guide counts disagree.</p>
            </article>
            <article className={styles.card}>
              <h3>2. Select guide_A</h3>
              <p>Its counts are 1, 0 and 1. Inspect how the matching rules account for the difference.</p>
            </article>
            <article className={styles.card}>
              <h3>3. Load read decisions</h3>
              <p>From the downloaded example, attach <code>bundle/read_changes.tsv</code>.
                Its identity is checked locally before any records appear.</p>
            </article>
          </div>
          <p className={styles.note}>The report stays on your machine. A matching checksum is not biological validation.</p>
        </section>
        <section className={styles.section}>
          <h2>Example results</h2>
          <p>
            This synthetic fixture contains close target sequences, a duplicate
            sequence with a distinct ID, substitutions, a literal N, an
            unmatched read and a short read. The values below are generated with
            DotMatch’s native Hamming matcher and checked in the test suite.
          </p>
          <div
            className={styles.scroll}
            tabIndex={0}
            role="region"
            aria-label="All policy outcomes"
          >
            <table className={styles.table}>
              <caption>All nine input records under each policy</caption>
              <thead>
                <tr>
                  <th scope="col">Policy</th>
                  <th scope="col">Unique</th>
                  <th scope="col">Ambiguous</th>
                  <th scope="col">Unmatched</th>
                  <th scope="col">Invalid</th>
                </tr>
              </thead>
              <tbody>
                {modes.map(([mode, label]) => (
                  <tr key={mode}>
                    <th scope="row">{label}</th>
                    {(["unique", "ambiguous", "none", "invalid"] as const).map(
                      (state) => (
                        <td key={state}>{demo.outcomes[mode][state]}</td>
                      ),
                    )}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <p className={styles.note}>
            Exact and radius-one each produce three unique assignments. But{" "}
            {demo.changed_reads} reads change outcome somewhere across the three
            policies. Equal totals do not establish equivalent count matrices—or
            equivalent read assignments.
          </p>
          <div
            className={styles.scroll}
            tabIndex={0}
            role="region"
            aria-label="Per-read policy decisions"
          >
            <table className={styles.table}>
              <caption>
                The underlying read decisions; no records omitted
              </caption>
              <thead>
                <tr>
                  <th scope="col">Read</th>
                  {modes.map(([mode, label]) => (
                    <th key={mode} scope="col">
                      {label}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {demo.records.map((record) => (
                  <tr key={record.id}>
                    <th scope="row">{record.id}</th>
                    {modes.map(([mode]) => (
                      <td key={mode}>
                        {record.calls[mode].target_id ??
                          record.calls[mode].status}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </section>
        <section className={styles.section}>
          <h2>Matching rules</h2>
          <div className={styles.grid}>
            <article className={styles.card}>
              <h3>Exact</h3>
              <p>
                Count a read only when exactly one target has the same sequence.
                Duplicate sequences under distinct IDs remain ambiguous.
              </p>
            </article>
            <article className={styles.card}>
              <h3>Radius one</h3>
              <p>
                Count only when one target is within one substitution. An exact
                read can become ambiguous when another target is a single
                substitution away.
              </p>
            </article>
            <article className={styles.card}>
              <h3>Best distance, k=1</h3>
              <p>
                Count when one target is nearest, allowing at most one
                substitution. Exact matches beat one-mismatch alternatives;
                equal-distance ties remain ambiguous.
              </p>
            </article>
          </div>
          <p className={styles.note}>
            None of these rules establishes a read’s true biological origin.
            This is a software and method-selection example, not an accuracy
            benchmark or an automatic recommendation to use a more permissive
            policy.
          </p>
        </section>
        <section className={styles.section}>
          <h2>Run on local data</h2>
          <p>
            <code>dotmatch sensitivity</code> is included in published {publishedVersion}.
            It produces the count tables and a static report. The interactive
            review on this page is a website example; it does not replace the
            local command's report.
          </p>
          <pre className={styles.code}>
            <code>{`python3 -m pip install dotmatch==${publishedVersion}\n\ndotmatch sensitivity \\\n  --targets guides.tsv \\\n  --reads sample.fastq.gz \\\n  --target-start 23 \\\n  --target-length 20 \\\n  --sample-label sample_1 \\\n  --out-dir sensitivity/`}</code>
          </pre>
          <p>
            Open <code>sensitivity/report.html</code>. The bundle contains three
            raw count tables, per-guide deltas, state transitions, sample QC,
            input and artifact checksums, and a machine-readable summary. Add{" "}
            <code>--write-read-changes</code> to record changed read IDs and
            calls without copying raw sequences.
          </p>
          <details>
            <summary>Open an existing run in the interactive viewer</summary>
            <p>
              From a current <a href={`${repoUrl}/tree/main`}>source checkout</a>,
              the standard-library renderer opens a completed v1 bundle without
              installing the native engine or rerunning your reads. The destination
              must be a new file; your original analysis is unchanged.
            </p>
            <pre className={styles.code}><code>{`python3 python/dotmatch/sensitivity_review.py --bundle sensitivity/ --out review.html`}</code></pre>
          </details>
          <p>
            The implementation reuses one native index and reads the FASTQ once.
            It compares the same fixed windows and does not change your baseline
            analysis, infer offsets or select an assignment policy.
          </p>
          <div className={styles.actions}>
            <a
              className={styles.primary}
              href={`${repoUrl}/tree/main/examples/assignment_sensitivity`}
            >
              Inspect the example inputs
            </a>
            <a className={styles.secondary} href={`${docsUrl}sensitivity.html`}>
              Command and output contract
            </a>
          </div>
          <p>
            Evaluating this on a real workflow? <a href={`${repoUrl}/issues/82`}>Share a de-identified result or reproducible discrepancy</a>.
            Do not post private reads, sample identifiers or unpublished libraries.
          </p>
        </section>
      </main>
      <ResearchFooter />
    </>
  );
}

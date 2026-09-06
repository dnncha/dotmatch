import demo from "../../public/assignment-demo.json";
import { ResearchHeader, ResearchFooter } from "../research-shell";
import {
  docsUrl,
  repoUrl,
  sitePath,
  pageMetadata,
  publishedVersion,
} from "../site-metadata";
import styles from "../research.module.css";
export const metadata = pageMetadata(
  "Why CRISPR guide counts change with matching rules | DotMatch",
  "Work through a native-checked example of exact, radius-one and best-distance assignment. Inspect changed reads, per-guide counts and reproducible outputs.",
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
          <p className={styles.eyebrow}>A worked scientific example</p>
          <h1 className={styles.title}>
            A matching percentage is not the whole result.
          </h1>
          <p className={styles.lede}>
            Exact and radius-one matching can assign the same number of reads
            while changing which guides receive the counts. Inspect the
            decisions, not just the total.
          </p>
        </header>
        <section className={styles.section}>
          <h2>Nine reads. Three explicit rules.</h2>
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
          <h2>Why do the policies disagree?</h2>
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
          <h2>Run the sensitivity workflow on your data.</h2>
          <p>
            The new <code>dotmatch sensitivity</code> command is in the source
            tree for the next release, not in published {publishedVersion}.
            Install from a reviewed source checkout before using it.
          </p>
          <pre className={styles.code}>
            <code>{`# From the DotMatch source checkout:\npython3 -m pip install .\n\ndotmatch sensitivity \\\n  --targets guides.tsv \\\n  --reads sample.fastq.gz \\\n  --target-start 23 \\\n  --target-length 20 \\\n  --sample-label sample_1 \\\n  --out-dir sensitivity/`}</code>
          </pre>
          <p>
            Open <code>sensitivity/report.html</code>. The bundle contains three
            raw count tables, per-guide deltas, state transitions, sample QC,
            input and artifact checksums, and a machine-readable summary. Add{" "}
            <code>--write-read-changes</code> to record changed read IDs and
            calls without copying raw sequences.
          </p>
          <p>
            The implementation reuses one native index and reads the FASTQ once.
            It compares the same fixed windows and does not change your baseline
            analysis, infer offsets or select a winning policy.
          </p>
          <div className={styles.actions}>
            <a
              className={styles.primary}
              href={`${repoUrl}/tree/main/examples/assignment_sensitivity`}
            >
              Reproduce this fixture
            </a>
            <a className={styles.secondary} href={`${docsUrl}sensitivity.html`}>
              Command and output contract
            </a>
          </div>
        </section>
      </main>
      <ResearchFooter />
    </>
  );
}

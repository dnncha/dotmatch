import { docsUrl, siteAsset } from "./site-metadata";
import styles from "./research.module.css";

export function GuideCounterBenchmark({ showGraphs = false }: { showGraphs?: boolean }) {
  return (
    <section id="guide-counter-benchmark" className={styles.section} aria-labelledby="guide-counter-title">
      <p className={styles.eyebrow}>Current source benchmark · One mismatch, no indels</p>
      <h2 id="guide-counter-title">DotMatch versus guide-counter</h2>
      <dl className={styles.stats}>
        <div className={styles.stat}><dt>Faster complete commands</dt><dd className={styles.benchmarkValue}>1.6–5.5×</dd></div>
        <div className={styles.stat}><dt>Lower peak memory</dt><dd className={styles.benchmarkValue}>11.3–14.4×</dd></div>
        <div className={styles.stat}><dt>Guide/sample count differences</dt><dd className={styles.benchmarkValue}>0</dd></div>
      </dl>
      <p>
        The benchmarked source checkout&apos;s <code>dotmatch guide-counter count</code> returned
        identical full count matrices in every paired comparison with unmodified guide-counter 0.1.3.
        The performance workloads use controlled 100k/1M-read FASTQs against the public
        87,437-guide Yusa library, with one or four samples.
      </p>
      <p className={styles.caption}>
        Five paired repeats, one CPU thread. Guide-counter is faster in the tested exact-mode
        million-read cases. Full experimental-screen performance remains unmeasured;
        these performance results do not establish biological accuracy.
      </p>
      <div className={styles.actions}>
        <a className={styles.secondary} href={`${docsUrl}benchmarks/guide_counter/README.html`}>
          Explore graphs, protocol and raw results
        </a>
      </div>
      {showGraphs && <>
        <figure className={styles.figure}>
          <a href={siteAsset("benchmarks/guide_counter_throughput.svg")}>
            <img src={siteAsset("benchmarks/guide_counter_throughput.svg")} width={1100} height={480}
              loading="lazy" alt="Complete-command throughput: DotMatch leads in one-mismatch workflows; guide-counter leads in exact-mode million-read cases." />
          </a>
          <figcaption>Complete commands include indexing, offset detection, gzip parsing, counting and all three output files. Bars show medians; whiskers show the observed five-run range.</figcaption>
        </figure>
        <figure className={styles.figure}>
          <a href={siteAsset("benchmarks/guide_counter_memory.svg")}>
            <img src={siteAsset("benchmarks/guide_counter_memory.svg")} width={1100} height={480}
              loading="lazy" alt="Peak memory: DotMatch uses about 37–47 MiB versus guide-counter’s 529 MiB in the tested one-mismatch workflows." />
          </a>
          <figcaption>Peak resident memory for the same commands and inputs. Both exact and one-mismatch modes are shown.</figcaption>
        </figure>
      </>}
    </section>
  );
}

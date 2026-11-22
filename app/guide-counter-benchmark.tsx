import { docsUrl, siteAsset } from "./site-metadata";
import styles from "./research.module.css";

export function GuideCounterBenchmark({
  showGraphs = false,
}: {
  showGraphs?: boolean;
}) {
  return (
    <section
      id="guide-counter-benchmark"
      className={styles.section}
      aria-labelledby="guide-counter-title"
    >
      <h2 id="guide-counter-title">
        Faster one-mismatch counting, identical counts
      </h2>
      <p>
        In paired benchmarks on the source checkout,{" "}
        <code>dotmatch guide-counter count</code> matched unmodified{" "}
        <code>guide-counter</code> 0.1.3 across full guide-by-sample count
        matrices on controlled 100k-read and 1M-read FASTQs against the
        87,437-guide Yusa library (one or four samples).
      </p>
      <dl className={styles.stats}>
        <div className={styles.stat}>
          <dt>Faster complete commands</dt>
          <dd className={styles.benchmarkValue}>1.6–5.5×</dd>
        </div>
        <div className={styles.stat}>
          <dt>Lower peak memory</dt>
          <dd className={styles.benchmarkValue}>11.3–14.4×</dd>
        </div>
        <div className={styles.stat}>
          <dt>Guide/sample count differences</dt>
          <dd className={styles.benchmarkValue}>0</dd>
        </div>
      </dl>
      <p className={styles.caption}>
        Five paired repeats on one CPU thread. <code>guide-counter</code> is
        faster in the tested exact-mode 1M-read cases; full experimental-screen
        performance and biological accuracy require separate validation.
      </p>
      <div className={styles.actions}>
        <a
          className={styles.secondary}
          href={`${docsUrl}benchmarks/guide_counter/README.html`}
        >
          View graphs, protocol, and raw CSVs
        </a>
      </div>
      {showGraphs && (
        <>
          <figure className={styles.figure}>
            <a href={siteAsset("benchmarks/guide_counter_throughput.svg")}>
              <img
                src={siteAsset("benchmarks/guide_counter_throughput.svg")}
                width={1100}
                height={480}
                loading="lazy"
                alt="Complete-command throughput: DotMatch leads in one-mismatch workflows; guide-counter leads in exact-mode 1M-read cases."
              />
            </a>
            <figcaption>
              Complete commands include indexing, offset detection, gzip
              parsing, counting, and all three output files. Bars show medians;
              whiskers show the observed five-run range.
            </figcaption>
          </figure>
          <figure className={styles.figure}>
            <a href={siteAsset("benchmarks/guide_counter_memory.svg")}>
              <img
                src={siteAsset("benchmarks/guide_counter_memory.svg")}
                width={1100}
                height={480}
                loading="lazy"
                alt="Peak memory: DotMatch uses about 37–47 MiB versus 529 MiB for guide-counter in the tested one-mismatch workflows."
              />
            </a>
            <figcaption>
              Peak resident memory for the same commands and inputs across both
              exact and one-mismatch modes.
            </figcaption>
          </figure>
        </>
      )}
    </section>
  );
}

"use client";
import { useState } from "react";
import demo from "../public/assignment-demo.json";
import styles from "./home.module.css";

const examples = [
  ["exact_isolated", "Exact"],
  ["one_mismatch", "One mismatch"],
  ["two_candidates", "Ambiguous"],
  ["unmatched", "Unmatched"],
  ["short", "Short read"],
] as const;

const explanations: Record<string, string> = {
  unique: "Matches one target within one substitution, adding one count to that target.",
  ambiguous:
    "Matches multiple targets within one substitution, so no count is added.",
  none: "Matches no target within one substitution; recorded as unmatched.",
  invalid:
    "Shorter than the required 20-base window; recorded as invalid.",
};

export function AssignmentDemo() {
  const [selected, setSelected] = useState<string>("one_mismatch");
  const record = demo.records.find((row) => row.id === selected)!;
  const result = record.calls.radius_k1;
  const target = demo.targets.find((row) => row.id === result.target_id);

  return (
    <div
      className={styles.instrument}
      aria-label="Interactive synthetic assignment example"
    >
      <div className={styles.instrumentHeading}>
        <span>Read matcher</span>
        <span className={styles.instrumentBadge}>Hamming k=1 (radius)</span>
      </div>
      <div
        className={styles.exampleControls}
        role="group"
        aria-label="Choose an example read"
      >
        {examples.map(([id, label]) => (
          <button
            key={id}
            type="button"
            aria-pressed={selected === id}
            onClick={() => setSelected(id)}
          >
            {label}
          </button>
        ))}
      </div>
      <div
        className={styles.sequenceArea}
        aria-live="polite"
        aria-atomic="true"
      >
        <div className={styles.sequenceBox}>
          <span className={styles.sequenceLabel}>Observed read window</span>
          <code className={styles.sequence}>
            {record.sequence.split("").map((base, i) => (
              <span
                key={i}
                className={
                  target && base !== target.sequence[i]
                    ? styles.mismatch
                    : undefined
                }
              >
                {base}
              </span>
            ))}
          </code>
        </div>
        <div className={styles.matchRule}>
          <span>
            {result.status === "invalid"
              ? "20-base window unavailable"
              : `${record.candidate_ids.length} compatible target${record.candidate_ids.length === 1 ? "" : "s"}`}
          </span>
          <span aria-hidden="true">↓</span>
        </div>
        <div className={styles.verdict}>
          <span className={styles.status} data-status={result.status}>
            {result.status === "none" ? "unmatched" : result.status}
          </span>
          <strong>
            {result.target_id ??
              (result.status === "ambiguous"
                ? "Not assigned"
                : "No count added")}
          </strong>
        </div>
        <p className={styles.explanation}>{explanations[result.status]}</p>
      </div>
      <p className={styles.instrumentFoot}>
        Synthetic fixture verified against the native matcher.{" "}
        <a href="#evidence">Inspect inputs and tests</a>
      </p>
    </div>
  );
}

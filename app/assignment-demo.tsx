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
  unique: "One target within one substitution. Adds one count to that target.",
  ambiguous:
    "Multiple targets within one substitution. No count is added.",
  none: "No target within one substitution. Recorded as unmatched.",
  invalid:
    "The read is shorter than the required 20-base window. Recorded as invalid.",
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
        <span>Read assignment example</span>
        <span>Radius-one Hamming</span>
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
        Synthetic example, checked against the native matcher.{" "}
        <a href="#evidence">Example inputs and tests</a>
      </p>
    </div>
  );
}

"use client";
import { useEffect, useRef, useState } from "react";
import { auditLibrary, parseLibrary, type SafetyReport } from "../../../lib/library-safety";
import styles from "../../research.module.css";
const example = "target_id\tsequence\nsynthetic_a\tACGTACGT\nsynthetic_b\tACGTTCGA\nsynthetic_c\tTGCATGCA";
export function LibraryExplorer() {
  const [text, setText] = useState("");
  const [report, setReport] = useState<SafetyReport | null>(null);
  const [error, setError] = useState("");
  const [status, setStatus] = useState("");
  const [busy, setBusy] = useState(false);
  const generation = useRef(0);
  useEffect(() => () => { generation.current++; }, []);
  function replaceText(value: string) {
    generation.current++; setText(value); setReport(null); setError(""); setStatus(""); setBusy(false);
  }
  function cancel() { generation.current++; setBusy(false); setStatus("Audit cancelled. No partial result is presented."); }
  async function run() {
    const id = ++generation.current;
    setReport(null); setError(""); setStatus("Checking the complete input…"); setBusy(true);
    try {
      const input = parseLibrary(text), iterator = auditLibrary(input);
      for (;;) {
        if (generation.current !== id) { iterator.return(undefined as never); return; }
        const next = iterator.next();
        if (next.done) { setReport(next.value); setStatus(`Complete: all ${next.value.target_count.toLocaleString()} targets checked.`); break; }
        setStatus(`${next.value.phase === "indexing" ? "Enumerating target neighbourhoods" : "Reviewing observations"}: ${next.value.completed.toLocaleString()} / ${next.value.total.toLocaleString()}`);
        await new Promise<void>(resolve => setTimeout(resolve, 0));
      }
    } catch (cause) {
      if (generation.current === id) { setError(cause instanceof Error ? cause.message : "The audit could not complete. No partial result is presented."); setStatus(""); }
    } finally { if (generation.current === id) setBusy(false); }
  }
  function download() {
    if (!report) return;
    const url = URL.createObjectURL(new Blob([JSON.stringify(report, null, 2) + "\n"], { type: "application/json" }));
    const anchor = document.createElement("a"); anchor.href = url; anchor.download = "dotmatch-library-audit.json";
    document.body.appendChild(anchor); anchor.click(); anchor.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  const rows = report ? [...report.targets].sort((a, b) => Number(b.exact_ambiguous_k1) - Number(a.exact_ambiguous_k1) || b.ambiguous_single_substitutions - a.ambiguous_single_substitutions || (a.id < b.id ? -1 : a.id > b.id ? 1 : 0)).slice(0, 100) : [];
  return <section className={styles.section} aria-labelledby="checker-title">
    <h2 id="checker-title">Check your target library</h2>
    <form className={styles.form} onSubmit={event => { event.preventDefault(); void run(); }}>
      <label htmlFor="library-input">Target sequences</label>
      <textarea id="library-input" className={styles.input} value={text} onChange={event => replaceText(event.target.value)} placeholder={example} spellCheck={false} autoCapitalize="off" autoCorrect="off" aria-describedby={error ? "library-help library-error" : "library-help"} aria-invalid={Boolean(error)} />
      <p id="library-help" className={styles.help}>Paste one sequence per line, or an unquoted two-column TSV/CSV: target_id and sequence. Optional header. Use 8–32 A/C/G/T bases, equal length, at most 2,000 targets. This tool makes no sequence-upload or storage requests.</p>
      <div className={styles.actions}><button type="submit" className={styles.primary} disabled={busy || !text.trim()}>{busy ? "Checking library…" : "Check complete library"}</button><button type="button" className={styles.secondary} onClick={() => replaceText(example)}>Load synthetic example</button>{busy ? <button type="button" className={styles.secondary} onClick={cancel}>Cancel</button> : <button type="button" className={styles.secondary} onClick={() => replaceText("")} disabled={!text}>Clear</button>}</div>
      {error && <p id="library-error" role="alert" className={styles.error}>{error}</p>}
      <p className={styles.status} role="status" aria-live="polite">{status}</p>
    </form>
    {report && <div aria-labelledby="result-title"><h2 id="result-title">Library geometry, not an error-rate estimate.</h2>
      <dl className={styles.stats}><div className={styles.stat}><dt>Targets checked</dt><dd>{report.target_count.toLocaleString()}</dd></div><div className={styles.stat}><dt>Targets with ambiguous exact reads at radius one</dt><dd>{report.targets_with_ambiguous_exact_reads.toLocaleString()}</dd></div><div className={styles.stat}><dt>Distinct ambiguous observations</dt><dd>{report.ambiguous_observations.toLocaleString()}</dd></div></dl>
      <p className={styles.readout}>{report.ambiguous_observations > 0 ? "Some target neighbourhoods overlap. Review the affected IDs and the matching policy before enabling correction." : "No overlapping radius-one neighbourhoods were found in the supplied library. This does not establish biological accuracy or cover errors outside this model."}</p>
      <p className={styles.help}>Hamming radius 1 · supplied orientation · radius ambiguity policy · {report.sequence_length} bases · {report.distinct_observations.toLocaleString()} distinct exact or one-substitution observations.</p>
      <div className={styles.actions}><button type="button" className={styles.secondary} onClick={download}>Save complete JSON report</button></div><p className={styles.caption}>The report includes your target IDs and sequences. Review it before sharing. No sequences are sent to DotMatch by this tool.</p>
      <div className={styles.scroll} tabIndex={0} role="region" aria-label="Per-target collision audit"><table className={styles.table}><caption>{report.target_count > 100 ? `Showing the first 100 of ${report.target_count} checked targets, ordered by ambiguity. All target results are in the JSON report.` : "All checked targets, ordered by ambiguity"}</caption><thead><tr><th scope="col">Target</th><th scope="col">Sequence</th><th scope="col">Exact read ambiguous at k=0</th><th scope="col">Exact read ambiguous at k=1</th><th scope="col">Ambiguous single substitutions</th></tr></thead><tbody>{rows.map(row => <tr key={row.id}><th scope="row">{row.id}</th><td><code>{row.sequence}</code></td><td>{row.exact_ambiguous_k0 ? "Yes" : "No"}</td><td>{row.exact_ambiguous_k1 ? "Yes" : "No"}</td><td>{row.ambiguous_single_substitutions} of {row.possible_single_substitutions}</td></tr>)}</tbody></table></div>
      {report.witnesses.length > 0 && <details className={styles.detail}><summary>Inspect ambiguous observation examples</summary><p>Each sequence below fits multiple target IDs under the stated radius rule. Up to 12 observation examples and 8 IDs per example are displayed; totals and per-target results are calculated over the complete input.</p><div className={styles.scroll} tabIndex={0} role="region" aria-label="Ambiguous observations"><table className={styles.table}><thead><tr><th scope="col">Observed sequence</th><th scope="col">Candidate targets</th><th scope="col">Target IDs</th></tr></thead><tbody>{report.witnesses.map(witness => <tr key={witness.observation}><td><code>{witness.observation}</code></td><td>{witness.candidate_count}</td><td>{witness.target_ids.join(", ")}{witness.ids_truncated ? " …" : ""}</td></tr>)}</tbody></table></div>{report.witnesses_truncated && <p>Examples are truncated; the aggregate audit is not.</p>}</details>}
    </div>}
  </section>;
}

/** Exhaustive Hamming-radius-one audit; not the native matcher or a biological error-rate model. */
export const LIMITS = Object.freeze({ targets: 2000, minLength: 8, maxLength: 32, characters: 400000 });
export type Target = { id: string; sequence: string };
export type Progress = { phase: "indexing" | "reviewing"; completed: number; total: number };
export type TargetAudit = Target & { exact_ambiguous_k0: boolean; exact_ambiguous_k1: boolean; ambiguous_single_substitutions: number; possible_single_substitutions: number };
export type SafetyReport = {
  schema_version: "dotmatch.library-safety.v1"; algorithm_version: "1.0.0"; evidence_level: "exact_combinatorial_audit";
  metric: "hamming"; radius: 1; ambiguity_policy: "radius"; orientation: "as_supplied";
  target_count: number; sequence_length: number; distinct_observations: number; ambiguous_observations: number;
  targets_with_ambiguous_exact_reads: number; targets_with_ambiguous_substitutions: number; targets: TargetAudit[];
  witnesses: { observation: string; candidate_count: number; target_ids: string[]; ids_truncated: boolean }[];
  witnesses_truncated: boolean; limitations: string[];
};
function validatedTargets(input: readonly Target[]): Target[] {
  if (!input.length) throw new Error("Add at least one target sequence.");
  if (input.length > LIMITS.targets) throw new Error(`This browser audit accepts at most ${LIMITS.targets} targets. Use dotmatch audit for a larger library; do not sample it to claim library-wide safety.`);
  const ids = new Set<string>();
  const length = input[0].sequence.length;
  const targets = input.map(({ id, sequence }) => {
    if (!id || id.length > 120 || !/^[\x20-\x7e]+$/.test(id)) throw new Error("Target IDs must contain 1–120 printable ASCII characters.");
    if (ids.has(id)) throw new Error(`Duplicate target ID: ${id}. Give each row a distinct ID; duplicate sequences may remain for review.`);
    ids.add(id);
    if (!/^[ACGT]+$/.test(sequence)) throw new Error(`${id}: use A, C, G and T only. N and other ambiguity codes are not silently discarded.`);
    if (sequence.length < LIMITS.minLength || sequence.length > LIMITS.maxLength) throw new Error(`${id}: sequences must be 8–32 bases long.`);
    if (sequence.length !== length) throw new Error("All sequences must have the same length for this Hamming audit. Use the native workflow for other matching rules.");
    return { id, sequence };
  });
  return targets.sort((a, b) => a.id < b.id ? -1 : a.id > b.id ? 1 : 0);
}
/** One sequence per line or unquoted two-column TSV/CSV, with an optional header. */
export function parseLibrary(text: string): Target[] {
  if (text.length > LIMITS.characters) throw new Error("Input is too large for this browser audit. Use a two-column library with at most 2,000 targets.");
  const lines = text.replace(/^\uFEFF/, "").split(/\r\n|\n|\r/).filter(line => line.trim().length > 0);
  if (!lines.length) throw new Error("Paste a library or load the synthetic example.");
  const targets: Target[] = [];
  for (const [index, line] of lines.entries()) {
    if (line.includes('"')) throw new Error("Quoted CSV is not supported here. Export an unquoted target_id / sequence TSV, or use the native workflow.");
    const fields = line.includes("\t") ? line.split("\t") : line.includes(",") ? line.split(",") : [line];
    const cells = fields.map(field => field.trim());
    if (index === 0 && ((cells.length === 1 && /^(sequence|seq)$/i.test(cells[0])) || (cells.length === 2 && /^(target_id|id|guide_id|barcode_id)$/i.test(cells[0]) && /^(sequence|seq)$/i.test(cells[1])))) continue;
    if (cells.length !== 1 && cells.length !== 2) throw new Error(`Row ${index + 1}: use one sequence or exactly two columns: target_id and sequence.`);
    targets.push({ id: cells.length === 2 ? cells[0] : `target_${String(targets.length + 1).padStart(4, "0")}`, sequence: cells[cells.length - 1].toUpperCase() });
    if (targets.length > LIMITS.targets) throw new Error(`This browser audit accepts at most ${LIMITS.targets} targets. Use dotmatch audit for the complete library.`);
  }
  return validatedTargets(targets);
}
/** Yields progress so browser callers can return control to the event loop. */
export function* auditLibrary(input: readonly Target[]): Generator<Progress, SafetyReport, unknown> {
  const targets = validatedTargets(input), length = targets[0].sequence.length;
  const owners = new Map<string, number[]>(), exactOwners = new Map<string, number>();
  const add = (sequence: string, owner: number) => { const previous = owners.get(sequence); if (previous) previous.push(owner); else owners.set(sequence, [owner]); };
  yield { phase: "indexing", completed: 0, total: targets.length };
  for (let i = 0; i < targets.length; i++) {
    const { sequence } = targets[i];
    exactOwners.set(sequence, (exactOwners.get(sequence) ?? 0) + 1); add(sequence, i);
    for (let position = 0; position < length; position++) for (const base of "ACGT") if (base !== sequence[position]) add(sequence.slice(0, position) + base + sequence.slice(position + 1), i);
    if ((i + 1) % 32 === 0) yield { phase: "indexing", completed: i + 1, total: targets.length };
  }
  const rows: TargetAudit[] = targets.map(target => ({ ...target, exact_ambiguous_k0: (exactOwners.get(target.sequence) ?? 0) > 1, exact_ambiguous_k1: (owners.get(target.sequence)?.length ?? 0) > 1, ambiguous_single_substitutions: 0, possible_single_substitutions: 3 * length }));
  let reviewed = 0, ambiguous = 0;
  const witnesses: SafetyReport["witnesses"] = [];
  yield { phase: "reviewing", completed: 0, total: owners.size };
  for (const [observation, candidates] of owners) {
    if (candidates.length > 1) {
      ambiguous++;
      for (const i of candidates) if (observation !== rows[i].sequence) rows[i].ambiguous_single_substitutions++;
      if (witnesses.length < 12) witnesses.push({ observation, candidate_count: candidates.length, target_ids: candidates.slice(0, 8).map(i => targets[i].id), ids_truncated: candidates.length > 8 });
    }
    if (++reviewed % 4096 === 0) yield { phase: "reviewing", completed: reviewed, total: owners.size };
  }
  return {
    schema_version: "dotmatch.library-safety.v1", algorithm_version: "1.0.0", evidence_level: "exact_combinatorial_audit", metric: "hamming", radius: 1, ambiguity_policy: "radius", orientation: "as_supplied",
    target_count: targets.length, sequence_length: length, distinct_observations: owners.size, ambiguous_observations: ambiguous,
    targets_with_ambiguous_exact_reads: rows.filter(row => row.exact_ambiguous_k1).length,
    targets_with_ambiguous_substitutions: rows.filter(row => row.ambiguous_single_substitutions > 0).length,
    targets: rows, witnesses, witnesses_truncated: ambiguous > witnesses.length,
    limitations: [
      "Library geometry only: these counts are not observed sequencing error rates, probabilities, FDR, or biological validation.",
      "Hamming radius 1, same-length ACGT sequences, supplied orientation, all supplied target rows; no indels, reverse-complement search, quality weighting, or best-distance policy.",
      "Duplicate sequences under different IDs remain separate candidates. No automatic deduplication or library repair is performed.",
      "A unique assignment can still be biologically wrong because of contamination, errors outside the model, or an incomplete target library.",
      "This standalone exhaustive implementation is not a native DotMatch benchmark or a certificate for an entire assay."
    ]
  };
}

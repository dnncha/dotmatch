# AR-001 amendments and execution log

## 2026-09-06 — implementation corrections

The first pilot run, 34026135700 at b6d4b2e982ff79911c9b451f1db1e3eff3532e6e, independently validated both 100,000-record Yusa evaluation prefixes and then stopped when the research parser treated Brunello's `Seq` column header as a sequence. Source inspection in run 34026360711 established that the unmodified Brunello source has header `sgRNAID, Seq, gene` (tab-delimited), exactly 77,441 data rows, and all sequences are valid 20-base ACGT strings. The parser now accepts the verified `Seq` header and includes a regression test. No target rows were dropped, trimmed, replaced or changed. This was a harness bug, not a reference-data defect or a DotMatch defect.

The first forensic run also stopped because guide-counter 0.1.3 does not implement `--version`. Version provenance now records that unsupported probe and validates the installed version through Cargo's installation registry, the installed binary SHA-256 and downloaded crate source hashes. Counting commands and scientific rules are unchanged.

Failed runs remain part of the audit trail. Initial Yusa outcomes were already inspected when these corrections were made. Neither correction was selected to improve the biological results. The original PROTOCOL.md is unchanged.

## 2026-09-06 — full-archive Yusa confirmation, extension AR-001-F

Frozen before inspection of any newly computed full-archive results. Motivation: nonrandom 100,000-read prefixes cannot establish full-run prevalence. Retrieve BOTH complete single-end archives ERR376998 and ERR376999, verify compressed bytes against ENA's advertised byte count and MD5, additionally record SHA-256, and require the actual parsed record count to equal archive metadata. Do not substitute prefixes if a full download fails.

The full-archive experiment remains a counting-methods study. It is NOT the biological inference stage B: these inputs provide one plasmid library and one cellular library, not independently replicated treatment groups. No gene-level p-values, FDR-controlled hit lists, biological-effect claims or false-assignment estimates will be generated from this extension.

Primary fixed-window comparison: exact, radius-one and best-distance Hamming at the already specified zero-based offset 23, length 19, same strand, full unchanged Yusa reference. Count all complete records, including invalid extractions. Independently reconstruct every guide count with a separate enumerated codeword table and validate that table against the query-neighbour oracle on exhaustive small fixtures and deterministic real windows. The three native DotMatch count matrices must agree exactly with the independent count vectors, including zeros and gene annotations. Required conservation checks: counts sum to unique calls and all read-state totals sum to records.

Separate multi-offset comparator lane: pinned guide-counter 0.1.3, Hamming-one mode, 100,000-record offset calibration, offset-min-fraction explicitly 0.0025. Reconstruct its offset selection and every guide count independently. Treat each selected-offset match as an event. Count distinct contributing read records separately; report zero/one/multiple events per read, same-guide repeats, different-guide repeats and same/cross-gene annotation cases. The original prefix lane remains distinct. Calibration reads are included in full-file event counting intentionally, matching the comparator; this is algorithmic reconstruction, not held-out biological validation.

Record full per-guide counts and aggregate event multiplicities. Save at most the first 1,000 multi-offset read ordinals as diagnostic traces, with the selection limit stated; do not represent those traces as a random sample. No raw reads are committed or published in artifacts. Full raw archives remain temporary working inputs. Results are complete only after all input checksums, record denominators, sample axes, annotations and independent count-equality gates pass. A failure must remain explicitly failed and preserve partial evidence.

Interpretation gate: an excess of matching events over read records establishes different counting units, not biological inaccuracy. Differences between a fixed-window tool and a multi-offset tool combine extraction coverage, match policy and event multiplicity. Do not call their entire count difference 'double counting' without reconciling these components. Library annotation equality does not establish true molecular origin.

# AR-001-B: Brunello historical count reconstruction

Frozen 6 September 2026 after viewing the six-sample fixed-window pilot and the successful Yusa prefix reconstruction, but before new Brunello comparator results. This is an exploratory implementation audit, not independent confirmatory evidence selected without prior observations.

Use the unchanged 77,441-guide corrected Brunello source and the first 100,000 records of each historical sample in PROTOCOL.md. These first-run prefixes include calibration records and are NOT full biological replicates or complete screen estimates. Preserve all configured source run metadata and each run's contribution count. Do not call the prefixes whole samples.

Reproduce the historical multi-offset command configuration: DotMatch 0.5.0 baseline, target start 20, length 20, Hamming k=1, best ambiguity policy, auto-offset range 20, calibration size 100,000, offset mode multi, offset minimum fraction 0.0025. Compare with guide-counter 0.1.3, calibration size 100,000 and fraction 0.0025. Identical library and per-sample FASTQs must feed both tools. If the original recorded aggregate (DotMatch 349,184; guide-counter 350,374; 255 differing guide totals) is not reproduced, record the actual result and investigate; do not tune parameters to force those numbers.

Independently reconstruct guide-counter's matching-event lookup and offset selection using the already tested forensic code. Validate its codeword lookup against query-enumeration calls on the first 1,000 canonical windows at start 23 in each prefix, including unmatched windows. Reconstruct every guide count in every sample, retain all per-sample discrepancies, and require zero reconstructed/comparator disagreements before marking reconstruction complete. A successful reconstruction says nothing about which original molecule generated a read.

Report matched read records separately from matching offset events, and classify multi-offset reads by same/different guide and same/cross-gene annotation. Counts remain separate by named sample axis. Post-processing may compare these events with DotMatch output and quantify source of differences; it may not infer biological accuracy from numerical equality.

No gene-level statistics or significance claims from these prefixes. A one-window lane is an extraction control, not a usable complete Brunello counting workflow when mixed guide offsets are present. Production code and historical artifacts remain untouched.

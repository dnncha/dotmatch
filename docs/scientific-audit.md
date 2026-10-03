# Scientific and implementation audit

Reviewed 2026-10-03. This audit covers the current working tree, including the
new counting comparison and CRISPRWorks Fit workflows. Performance superiority
on the tested counting workloads does not establish scientific SOTA for
gene-level inference.

## Fixed in this audit

| Issue | Fix and verification |
| --- | --- |
| MLE could overwrite a count/design/control input with its running manifest or another output. | Resolve existing input paths and reject final-output collisions before writing. Regressions cover all four advertised output suffixes and preserve original bytes. |
| Joint and calibration used predictable `.tmp` paths, which could consume an input or another run's staging file. | Reserve exclusive, unique sibling staging files. Calibration regression uses valid fit details named like the previous temporary output and verifies both input bytes and file provenance survive. |
| Path checks missed hard links to inputs. | Compare file identity as well as resolved paths; MLE hard-link and shared symlink/hard-link regressions preserve inputs. |
| Publication could replace a table before manifest staging failed. | Stage both files first, remove an old completion manifest before replacing the table, and publish the new manifest last. Fault-injection tests preserve the old bundle on staging failure and leave no misleading manifest after publication failure. |
| Successful joint fits were excluded by the MLE permutation guide-count cap. | Apply the cap only to MLE records, report no joint cap, and test 40- and 80-guide joint records. |
| Joint details without explicit model metadata were labeled as MLE effects with the wrong units. | Infer the default model and units from the validated effect schema; verify joint metadata in the CLI regression. |
| Fit CI used unittest discovery, which skipped the new pytest-style publication tests. | Install the package's test extra and run the full suite with pytest in the existing platform/Python matrix. |
| Joint output did not distinguish meeting the stopping rule from iteration exhaustion. | Record per-gene iterations, termination reason and final stopping-statistic change; expose validated CLI iteration/tolerance settings. Defaults retain upstream numerical parity. |
| Calibration accepted fits explicitly marked nonconverged. | Reject nonconverged or malformed diagnostics before publishing calibration results; legacy records without diagnostics remain supported with that limitation. |
| Calibration ignored the companion summary's recorded hash. | Verify recorded companion outputs, preserve them as inputs, and reject missing/tampered/malformed records and output collisions. This verifies bundle consistency, not biological validity or the authenticity of fit-detail contents. |
| A convergence flag alone could bypass calibration's diagnostics check. | Require positive iteration counts, finite nonnegative stopping changes and a matching termination reason; when declared stopping settings are present, reject contradictions and incomplete diagnostics. |
| The benchmark validation file could omit every figure while still passing. | Require the complete six-figure set, valid digest syntax, the recorded command count and zero-error summary; reject missing and unexpected artifacts. |
| Website and documentation headlines were not checked against raw ratios. | Derive expected rounded intervals from paired measurements and require them in all twelve known callout sources; a stale website headline regression exercises rejection. This is a presence check, not a general validator of every numeric statement. |
| Calibration could silently report no calls when control resolution made discoveries impossible. | Report optimistic directional/two-sided q-value bounds for the declared family and current eligible strata; warn when no hypothesis can reach the chosen threshold. A 20-control/500-extreme-test regression verifies the BY limitation against SciPy. P/q outputs are unchanged. |
| Native offset forwarding rounded thresholds and changed selected windows at a boundary. | Preserve 17 significant digits in both forwarding paths and QC metadata. Independent byte-oracle tests cover immediately below/equal/above 0.5 in exact and one-mismatch modes with plain/gzip reads. Rerun the complete paired benchmark for the corrected native source. |
| Benchmark statistics comparison accepted NaN differences and extended parsing overwrote duplicate guides. | Reject nonfinite statistics, incomplete/duplicate columns, duplicate extended rows and missing sample coverage. Add harness regressions and apply these checks in the refreshed paired run. |
| Cooperating Fit publishers could interleave a table and manifest at the same prefix. | Guard publication with a nonblocking OS lock; reject overlapping writers before replacement, retain the lock inode and rely on OS release after process exit. A separate-process test kills the holder and verifies subsequent publication succeeds. |
| Input provenance was recorded after fitting, allowing changed inputs to be mislabeled. | Parse joint count/map/control and calibration detail/control inputs from immutable byte payloads whose exact hashes are recorded. Recheck file hashes and output aliases before publication; changing-input and late-collision tests preserve earlier results. |
| MLE wrote tables directly and could overwrite an earlier completed run during a failed rerun. | Fit from private input snapshots into a same-directory staging area; require all outputs before locked publication, check original input hashes, and publish the completion manifest last. Fitting-failure and changed-input regressions preserve the prior bundle; mid-publication failure leaves no completion manifest. |
| Rerunning MLE without optional details left stale fit details from the earlier run. | Remove stale optional details only during successful publication, protecting that destination against input aliases even when detail output is disabled. Preserve CSV snapshot suffixes because upstream parsing depends on them. |
| Upstream MLE loading accepted nonfinite/negative counts and duplicate guide IDs. | Validate the immutable count snapshot before fitting: distinct nonempty columns, consistent rows, nonempty IDs, globally unique guide IDs and finite nonnegative counts. Regressions verify invalid counts never reach the fitter and prior outputs survive. |
| Benchmark evidence validation disappeared under `python -O`. | Replace assertions with explicit rejection branches. A subprocess regression proves optimized Python rejects a corrupted count comparison even when its CSV hash is updated. |

Input-preservation tests are in
`packages/crisprworks-fit/tests/test_publication.py`; benchmark validation tests
are in `python/tests/test_report_guide_counter.py`. Benchmark timings were not
changed by these fixes.

## Remaining priorities

1. **Filesystem crash durability.** MLE, joint and calibration now stage their
   required outputs and guard publication with OS locks. Failed fitting or
   staging preserves the prior completed bundle. Publication invalidates the
   old manifest before replacing files and installs a new manifest last; a
   failure partway through leaves an incomplete bundle, not a rolled-back one.
   This does not provide filesystem/machine crash durability or constrain
   unrelated programs writing the same paths. Consumers must verify completion
   manifests and output hashes; evaluate versioned run directories and durable
   fsync/rename ordering if crash recovery is required.
2. **Calibration power and validity.** The new preflight reports attainable
   resolution, but does not improve it. The current frontier panel still has
   no whole-gene BY calls. Evaluate larger independent control panels and
   end-to-end count-generating simulations for both power and error control.
   Switching to BH requires evidence for its assumptions; optimistic bounds
   are neither expected power nor an error-control guarantee.
3. **Explicit MLE fit status.** The joint exclusion bug is fixed. MLE still
   infers `skipped_fit` from its configured permutation guide-count cap rather
   than storing an explicit per-gene fit/permutation status. Carry those
   diagnostics in the fit details.
4. **Stopping-rule adequacy and legacy diagnostics.** New joint results expose
   numerical termination and calibration rejects explicit nonconvergence.
   Meeting the inherited bound-change threshold does not establish posterior
   accuracy. Evaluate stability under stricter settings and difficult count
   scenarios; older records without diagnostics cannot establish convergence.
5. **Independent biological validation.** Current gene-ranking evidence uses
   control-enriched subsets of two libraries and three cell lines. Extend to
   genome-wide independent screens, variable guide counts, inactive guides and
   copy-number confounding before claiming broad scientific superiority.
6. **Experimental counting performance.** The guide-counter speed headline is
   from controlled simulated reads. The cached experimental fixtures contain
   only 25 reads each. Full experimental FASTQs, quality/error strata and
   multiple thread counts are still needed. Exact-only million-read timings
   favor guide-counter and must remain visible.
7. **Claim synchronization beyond headline presence.** The gate now checks
   expected speed/memory intervals in all twelve known callout sources. Those
   literals still need a shared generated source. The presence check can miss
   an incorrect duplicate statement when a correct interval remains elsewhere
   in the same file; it does not validate graph alt text or every numeric claim.
8. **Metric-aware generic index construction.** The generic matching index
   still builds deletion storage and Hamming seeds without a caller-selected
   metric/radius. Inspect and benchmark immutable metric-specific construction
   to avoid unused allocation for exact-only and Hamming-only callers. The
   dedicated guide-counter path uses its own compact lookup, so its memory
   results do not establish generic-index memory efficiency.


The [counting report](benchmarks/guide_counter/README.md) gives the measured
performance scope and exact-mode losses. Fit biological evaluation and
calibration limits are documented in the package's `benchmarks/FRONTIER.md`
and `benchmarks/NULL-CALIBRATION.md`. These remaining issues are not represented
as resolved by the counting compatibility tests or speed graphs.

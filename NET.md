# Feedback net

Inner loop: `make inner` (about 15 s warm on the Linux worker; also the pre-commit hook)
Gate: `make gate` (84 s warm, about 8 min cold on the Linux worker; also the pre-push hook)
Scheduled: `make nightly` (weekly T3 scheduled task on a Linux worker, red files a GitHub issue)

Enable the hooks once per checkout with `make hooks`. GitHub Actions are billing-blocked and are not part of this net; local gate results with the commit SHA are the CI record.

## Signals by tier

| Tier | Signal | Command | Runs | Time | Status |
|---|---|---|---|---|---|
| edit | C banned-API and output-close rules, Makefile pipefail rule | `make lint-c` (`tools/lint_c.py`) | hook, inner, gate | 0.3 s | honest |
| edit | Strict C build: `-Wall -Wextra -Wpedantic -Wformat=2 -Wconversion -Wvla -Wimplicit-fallthrough -Wshadow -Werror` | `make strict` | hook, inner, gate | 40 s cold, 0 s warm | honest |
| edit | Python lint (typed sources full set, pyflakes everywhere) | `make lint` | hook, inner, gate | 1 s | honest |
| inner | C unit tests | `make test` | hook, inner, gate | 3 s warm | honest |
| inner | qdalign differential fuzz vs O(nm) reference, ASan/UBSan, 3k cases | `make fuzz-smoke` | hook, inner, gate | 3 s | honest |
| inner | FASTQ/CLI differential fuzz (plain vs gzip, count vs fastq-assign, lookup vs generic engine, lookup vs scan), ASan/UBSan, 80 seeds | `make fuzz-smoke` | hook, inner, gate | 2 s | honest |
| gate | NET.md wiring and rung check | `make netcheck` (`tools/netcheck.py`) | hook, gate | 0.1 s | honest |
| gate | CLI behaviour tests | `make cli-test` | hook, gate | 95 s | honest |
| gate | ASan/UBSan unit tests plus the CLI suites against the sanitizer binary | `make asan` | hook, gate | 3 min | honest |
| gate | mypy strict and pyright | `make typecheck` | hook, gate | 80 s | honest |
| gate | Python suite (2,160+ tests incl. readiness and packaging checks) | `make python-test` | hook, gate | 50-145 s | honest |
| release | Full pre-tag gate | `make pretag-ready` | release | 30+ min | honest |
| weekly | Long fuzz (400k qdalign cases, 6,000 FASTQ seeds, time-based seed) | `make fuzz-long` via `make nightly` | scheduled: T3 schedule_task on a Linux worker, Mondays 03:00 UTC, red files a GitHub issue | 25 min | honest |
| weekly | Mutation sensitivity of the scientific invariants | `make mutation-scientific` via `make nightly` | scheduled: T3 schedule_task on a Linux worker, Mondays 03:00 UTC, red files a GitHub issue | 5 min | honest |
| human | Distribution channel verification | `make distribution-channels` | human, per release | 1 min | red until the 0.8.0 PyPI wheel, GHCR tag and Bioconda recipe land |

## Risk surfaces

| Surface | Paths | Risk (L×I) | Covered by | Status |
|---|---|---|---|---|
| FASTQ reader and record validation | `src/qda.c` (`fastq_*`), `python/dotmatch/fastq_io.py` | high × high | one reader, fuzz-smoke, `tests/test_native_fastq_validation.py` | covered |
| Assignment engines (lookup, generic, scan, status shortcut) | `src/qda.c` (`assign_*`, `*_lookup_*`), `src/qdalign.c` | high × high | fuzz-smoke differential (both harnesses), `make test` | covered |
| Target, sample and barcode tables | `src/target_table.h`, `src/qda.c` (`read_table`, `read_samples_file`), `python/dotmatch/target_io.py` | med × high | lint-c fgets ban, cli-test, `python/tests/test_target_io.py` | covered |
| Output files and exit codes | `src/qda.c` (`open_output_file`, `close_output_file`) | med × high | lint-c close rule, /dev/full checks in history | covered |
| Allocation sizes | `src/qda.c`, `src/qdalign.c` | low × high | lint-c multiplied-size rule | covered |
| Quality-gated corrections | `src/qda.c` (`quality_allows_unique_correction`) | med × med | cli-test regression | late: regression test only |
| Paired reads | `src/qda.c` (`run_pair_count`), `python/dotmatch/assayruntime.py` | med × med | cli-test, `python/tests/test_assayruntime.py` | late: regression tests only |
| Python CLI fallback vs native CLI | `python/dotmatch/cli.py`, `python/dotmatch/core.py` | med × med | `python/tests/test_python_native_differential.py` (30 seeds across 14 metric/k/policy combinations) in `make python-test` | covered |
| Release records and evidence claims | `docs/distribution-release.json`, `CITATION.cff`, `README.md` | med × med | python-test readiness checks, pre-push gate | covered |

## Holes, ranked

1. No mutation testing of the C core (`mutation-scientific` mutates Python only). Next step: a small mutant runner over `src/qdalign.c` comparison operators scored by `make test fuzz-smoke`, weekly.
2. Coverage-guided fuzzing is unavailable on the worker toolchain (clang 14 without compiler-rt). The random harnesses stand in. Next step: a libFuzzer target over `fastq_next_record` on a machine with compiler-rt.
3. Python and native FASTQ/target parsers still differ on edge bytes: bare `\r` line endings (Python accepts), `\x1c` in headers, NBSP in IDs, reads over 8,191 bases (native rejects with a positioned error, Python accepts). Next step: decide one contract per case and add each to `tests/test_native_fastq_validation.py`.
4. `COVERAGE_MIN ?= 75` in the Makefile fails `make coverage` on a fitted number with no stated user-facing reason. It is only reached from `pretag-ready`; next step: print the number instead of failing.

## Escape log

| Date | Symptom | Root cause | Class | Caught by | Fixed (+siblings) | Class guard | Rung | Time to red | Proved |
|---|---|---|---|---|---|---|---|---|---|
| 2026-11-23 | main red: python-test (README evidence boundary, 0.8.0 DOI record) and pyright | copy and release edits pushed with no local gate run; Actions billing-blocked | unrun gate before push | dragnet recon | f43f3d7 | `.githooks/pre-push` runs `make gate`, `make hooks` enables it | 3 rule (hook) | 8 min | yes: gate red on 64e9c22 |
| 2026-11-23 | NUL byte silently truncated a plain FASTQ read; long lines split or passed depending on command | three line readers (fgets+strlen, getline, gzread) with different semantics | divergent duplicate readers | new fuzz harness (seed 11, 220) | 28889ac (sequence-only reader too) | plain FILE* path and duplicate record reader deleted; every input goes through gzopen and `fastq_next_record` | 2 deleted path | n/a | yes: fuzz red on f43f3d7 |
| 2026-11-23 | assign/match and crispr-count --samples made phantom reads and samples from long lines | fgets into fixed buffers with no newline check | fixed-buffer line reads | parser audit | d520d77 (BCL sample sheet, HTML preview) | `tools/lint_c.py` bans fgets/gets/sprintf/strcpy/strcat/strtok/ato* | 3 rule | 0.3 s | yes: planted fgets red |
| 2026-11-23 | count/fastq-assign/pair-count/demux/audit exited 0 on /dev/full | unchecked fclose on outputs, stdout never checked | dropped write errors | parser audit | f98c267 | lint rule: no fclose on an open_output_file handle, every handle reaches close_output_file; main checks stdout | 3 rule | 0.3 s | yes: planted fclose red; /dev/full exit 0 before, 1 after |
| 2026-11-23 | heap overflow with huge --auto-offset-sample | count * size wrapped in realloc | unchecked allocation product | parser audit | aee65ca (62 sites) | malloc_array/realloc_array; lint rule rejects multiplied malloc/realloc sizes | 3 rule | 0.3 s | yes: planted malloc(n * size) red |
| 2026-11-23 | status_stats AMBIGUOUS with wrong best distance; ascii lookup missed lowercase targets; empty-index SEGV | early exit before the remaining candidates' floor; one-sided case fold | engine shortcut disagrees with reference | qdalign differential harness | 70a4900, 198d21f | `tests/fuzz_qdalign_diff.c` in fuzz-smoke | 6 property or fuzz; above: the invariant is search-order correctness over the candidate space, which no type or lint rule can state | 11 s | yes: red on 28889ac and on a planted unconditional stop |
| 2026-11-23 | count lookup engine dropped N windows; indel ties classified differently per engine; 64-bit shift UB | lookup returned NONE instead of falling back; generic loop visited target_len-1 first | engine shortcut disagrees with reference | parser audit, then fuzz engine check | bf61100 | fuzz-smoke compares lookup engine with generic engine (forced by --assignments) | 6 property or fuzz; above: equivalence of two engines is behavioural, not a property of any one call site | 1 s | yes: red at seed 2 on a planted missing fallback |
| 2026-11-23 | Levenshtein --max-correction-qual rescued Q40 reads as deletions | deletion branch returned 1 with no quality check | quality gate bypass | parser audit | 27c08c1 | regression in `tests/test_cli_fastq.sh` | 7 regression test; above: which base quality justifies an indel is a domain rule with no type or lint encoding | 95 s | yes: red on f43f3d7 |
| 2026-11-23 | swapped, repeated or identical R1/R2 accepted | /1 and /2 stripped from either mate | mate identity unchecked | parser audit | 052dcc3 (AssayCode runtime sibling) | regressions in `tests/test_cli_fastq.sh`, `python/tests/test_assayruntime.py` | 7 regression test; above: mate order is a property of record pairs read at run time | 95 s | yes: red on 27c08c1 |
| 2026-11-23 | demux --barcode-length auto unique although a second barcode was inside the radius | per-length merge dropped earlier hits | lossy result merge | qdalign audit | c988e46 | regression in `tests/test_cli_fastq.sh` | 7 regression test; above: one remaining caller of the lossy merge is correct for the best policy, so the function cannot be deleted | 95 s | yes: unique before, ambiguous after |
| 2026-11-23 | a header row became a target (half headers natively and in Python; barcode headers in fastq-assign) | header detection required both columns; read_table never detected headers | header row read as data | qdalign and Python audits | 67134e9 | regressions in `tests/test_cli_fastq.sh`, `python/tests/test_target_io.py` | 7 regression test; above: header recognition is a vocabulary check on file content | 95 s | yes: red on f43f3d7 |
| 2026-11-23 | Python ambiguity flags, status-path radius, posterior with mixed lengths, readers raising EOFError | one-sided neighbour search, early-stop counts used for radius, whole library passed to a same-length model, missing exception mapping | Python helper semantics | Python audit | 21b3469, df9af2b | regressions in `python/tests/test_cli.py`, `test_matcher.py`, `test_io_boundaries.py` | 7 regression test; above: each is a single-call-site semantic rule; the parser-level gap is hole 1 | 50 s | yes: red with the fixes stashed |
| 2026-11-23 | 70a4900 left make test red; a piped fuzz check could not fail the gate | test output filtered through grep while iterating; `cmd \| tail -1` in a recipe | hidden check exit code | make inner | 198d21f, this commit | `.SHELLFLAGS := -o pipefail -ec` and a lint rule requiring it | 3 rule | 0.3 s | yes: lint red with pipefail removed; fuzz-smoke red on a planted stop |
| 2026-11-23 | pyflakes autofix removed `guide_counter_offset_args`, which a test reads as a module attribute of the benchmark script | `ruff --fix` cannot see attribute access through a dynamically loaded script | lost re-export | pre-push `make gate` (python-test) | this commit | explicit `name as name` re-export, which ruff treats as used; the gate test stays | 7 regression test; above: the use is a getattr on a module loaded by path, invisible to static analysis | 113 s | yes: gate red on cecfc26 |

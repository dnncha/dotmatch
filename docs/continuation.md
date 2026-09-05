# Next working session

This is an implemented research-alpha package, not a request to generate a new
scaffold. Read `BUILD_STATUS.md`, `AGENTS.md` and `roadmap.json` first.

Start with EW-001. In a local checkout with an authenticated GitHub CLI, run the
source inventory check and the tests. Create the standalone repository with
`python scripts/publish_github.py --public` only when `dnncha/editwitness` does not
already exist. Never merge the staging branch into DotMatch. If the target exists,
inspect it and continue its existing history rather than creating a replacement.

All six Linux/macOS/Windows jobs, strict mypy, schemas, distribution builds and
coverage passed in GitHub run 33998897914 for commit
`377c3bfafeb722bce3eac91d339dfec9d1028e68`; see `BUILD_STATUS.md`. The schema typing
and Windows UTF-8 test fixes are already included. Repeat the checks in the new
standalone repository after publication. Ruff was not run; do not claim otherwise.
Fix demonstrated failures with narrow patches and regression tests, not suppressed
checks. Record the actual commit and run links in reviewed status updates.

Confirm ownership/availability of the PyPI name before publishing. Configure a
reviewed trusted-publishing path and release only an explicitly alpha version.
There is no stored API key or assumed package-index account. Run the installed
wheel's demo and replay commands from outside the repository after building.

Then advance a bounded, unblocked roadmap item. Prioritize independent scientific
review, exact primer rematching, useful hypothesis generation and adjudicated
biological examples over more UI surfaces. External review and biological data
must be obtained, not fabricated. Keep the original model replayable, and bump the
model identifier when its semantics change.

Finish each iteration by reporting implemented changes, exact tests actually run,
remaining scientific limits and one highest-value next action. Update roadmap
statuses only with evidence. Do not turn software correctness into claims of
measured assay sensitivity or clone safety.

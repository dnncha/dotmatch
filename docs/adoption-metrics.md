# Adoption Metrics

DotMatch reports distribution health separately from evidence of real workflow
use. The machine-readable contract is [`adoption-metrics.json`](adoption-metrics.json)
and is checked by `make adoption-metrics-ready`.

## What the numbers mean

- Anaconda and PyPI downloads are package-distribution signals. They include
  mirrors, CI, repeated installs, and automated jobs; they do not identify
  active teams or successful runs.
- The north-star metric is completed independent evaluations using the intake,
  run, output review, and scorecard in [`pilot-program.md`](pilot-program.md).
- Retention is measured by repeat workflows, not by a second download.
- Ecosystem progress is an accepted external integration recorded in
  [`workflow-adoption.json`](workflow-adoption.json), with public evidence URLs.

## Monthly review

1. Record the Anaconda and PyPI distribution counts and the date collected.
2. Count completed evaluation scorecards, keeping private assay data out of the
   repository.
3. Count repeat workflows and accepted external integrations from their source
   records.
4. Compare the results with the rolling 90-day targets in the JSON contract.
5. Update positioning or onboarding only when the completed-evaluation pattern
   identifies a concrete friction point.

The project should not publish a unique-user, active-project, production-use,
or market-share claim from package download counts alone.

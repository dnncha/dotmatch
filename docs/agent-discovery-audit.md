# Agent discovery audit

This audit measures whether a coding or scientific agent can discover DotMatch,
select the correct entry point, identify its limitations, and verify an
installed workflow. It does not score scientific quality, adoption, popularity,
or unique users.

## Reproduce the measurements

Candidate worktree:

```bash
python3 scripts/check_agent_discovery.py --measure --json
```

Verified baseline commit `e0bedb8fb2bf514393c080c653960f17b6ae7856`:

```bash
git archive --format=tar --output=.agent-discovery-baseline.tar \
  e0bedb8fb2bf514393c080c653960f17b6ae7856
mkdir .agent-discovery-baseline
tar -xf .agent-discovery-baseline.tar -C .agent-discovery-baseline
python3 scripts/check_agent_discovery.py \
  --root .agent-discovery-baseline \
  --measure \
  --json
```

Currently deployed public surfaces:

```bash
python3 scripts/check_agent_discovery.py --live --json
```

The live mode sends read-only requests to the GitHub repository API, PyPI JSON
API, Read the Docs, and GitHub Pages. Network checks are deliberately separate
from the CI gate so provider outages do not make source validation
nondeterministic.

## Baseline and candidate

Measured on 2026-08-27:

| Surface | Baseline | Candidate | Meaning |
| --- | ---: | ---: | --- |
| Local agent-native discovery checks | 0/12 | 12/12 | The baseline lacked the newly defined task table, installed JSON route, capability schema, `llms.txt` files, synced public copies, structured agent links, PyPI agent URL, and clean-install workflow hook. |
| Currently deployed public checks | 4/7 | Not deployed | GitHub search metadata, PyPI scope metadata, Read the Docs onboarding, and GitHub Pages structured scope pass. Pages `llms.txt`, Read the Docs `llms.txt`, and the public capability manifest return 404 before merge/deployment. |

The 0/12 local baseline is not a claim that the earlier documentation was
empty. The 4/7 live score records its real strengths. The local rubric measures
the additional agent-native contracts introduced by this change and gives each
contract one equal, binary check.

## Local rubric

1. README task routing and Agent guide link.
2. Installed help task routing.
3. Installed `dotmatch capabilities --json` interface.
4. Canonical capability manifest.
5. Valid required intents, inputs, outputs, limitations, and evidence paths.
6. JSON Schema draft 2020-12 contract.
7. Concise `llms.txt`.
8. Self-contained `llms-full.txt`.
9. Byte-identical GitHub Pages, Read the Docs, and installed-package copies.
10. Web `describedby` link and `SoftwareApplication` feature list.
11. PyPI Agent guide project URL in source metadata.
12. Fresh-venv wheel and source-distribution FASTQ workflow gate.

`make agent-discovery-ready` fails when any final contract is absent, invalid,
or stale. `make repository-ready` calls the same gate.

## Clean-environment execution

`make python-package-test` passed on macOS on 2026-08-27. The gate built both
`dotmatch-0.2.2.tar.gz` and the universal macOS wheel, installed each artifact
into a separate new virtual environment, and verified:

- the package and native command report version `0.2.2`;
- `dotmatch capabilities --json` contains all six required task ids;
- a four-read, four-target `k=0` FASTQ count reports four unique assignments;
- the expected target row is present in the count table;
- the existing AssaySpec check, run, inference, autopsy, and CRISPR QC package
  smoke paths still complete.

CI repeats the package gate on Linux and macOS. This local result is package
execution evidence, not evidence that the unreleased command is already
available from PyPI.

## Intent coverage

The manifest routes seven representative intents. Every row includes queries,
an exact command template, required inputs, outputs, limitations,
documentation, and repository evidence.

| Manifest id | Entry point | Limitation that prevents over-routing |
| --- | --- | --- |
| `crispr-guide-counting` | `dotmatch crispr-count` | No downstream screen statistics |
| `inline-barcode-demultiplexing` | `dotmatch demux` | FASTQ input only; no basecalling |
| `feature-barcode-assignment` | `dotmatch count` | Per-read only; no cell/UMI or Cell Ranger quantification |
| `perturb-seq-guide-capture` | `dotmatch count` | Public evidence is single-guide extraction; no guide-per-cell or perturbation effects |
| `barcode-panel-design` | `dotmatch panel design` | Short barcode sets, not probe or full assay design |
| `known-target-fastq-matching` | `dotmatch count` | Finite known targets and a reviewed fixed window |
| `barcode-run-diagnosis` | `dotmatch barcode autopsy` | Suggestions still require assay-context review |

## Publication gate

The candidate does not change repository visibility or claim that its source
state is already live. After merge, the GitHub Pages and Read the Docs builds
must publish `llms.txt`, `llms-full.txt`, and the capability JSON at their
documented URLs. PyPI will expose the new Agent guide project URL and the
installed `capabilities` command only after a versioned package release passes
the existing distribution gates.

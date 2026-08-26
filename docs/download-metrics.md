# Download metrics

This snapshot records provider-reported package retrieval events across public
distribution channels. It does not claim that those events represent distinct
users. Distribution counts remain separate from discovery signals and evidence
of scientific use.

The current snapshot schema is version 2: it includes measurable GitHub
Release asset downloads and the observed BioContainers/Quay pull series
and Zenodo release-record downloads alongside Anaconda and PyPI events.

## Refresh the snapshot

From the repository root:

```bash
make download-metrics
```

This updates [`download-metrics.json`](download-metrics.json) from public
provider endpoints. The snapshot records the capture time, source URLs, channel
counts, package version/platform/Python breakdowns, the aggregate total,
and the interpretation boundary.

The repository also refreshes this snapshot weekly through the
`download-metrics` GitHub Actions workflow. A manual `make download-metrics`
run remains available for release and investigation checkpoints.

## Sources and boundaries

| Source | What is measured | What it cannot establish |
| --- | --- | --- |
| [Bioconda package API](https://api.anaconda.org/package/bioconda/dotmatch) | Cumulative downloads of published package artifacts | Distinct users, organizations, location, or runtime use |
| [Anaconda package metadata](https://api.anaconda.org/repocore/channels/bioconda/artifacts/conda/dotmatch) | Current package total and release metadata | Whether a download came from a person, CI, or a cache miss |
| [PyPI Stats](https://pypistats.org/packages/dotmatch) | Observed retained daily package downloads excluding mirrors, the current provider window, and recent monthly downloads | A provider lifetime total, distinct users, or the purpose of each install |
| [GitHub Releases](https://api.github.com/repos/dnncha/dotmatch/releases?per_page=100) | Public release-asset downloads by tag, asset, and broad platform | Distinct users, source-code use, or runtime use |
| [GitHub repository](https://github.com/dnncha/dotmatch) | Stars and forks as discovery context | Package downloads or scientific adoption |
| [GHCR container](https://github.com/dnncha/dotmatch/pkgs/container/dotmatch) | Public image availability | Pull totals are not exposed by the unauthenticated registry endpoint |
| [BioContainers / Quay](https://quay.io/repository/biocontainers/dotmatch) | Provider-reported image pulls in the current Quay stats window | Not a lifetime total; Quay does not attribute pulls by version, platform, client, or runtime use |
| [Zenodo release](https://doi.org/10.5281/zenodo.21511337) | Cumulative downloads of the published release record | Distinct users, client identity, or runtime use |

Anaconda’s artifact count is broken down by version, platform, and Python
build. GitHub Release assets are broken down by tag, asset, and broad platform.
Multiple artifacts and build-number updates make automated or matrix traffic
plausible, so the report labels these sources `mixed_unknown` rather than
converting events into user estimates.

The aggregate number sums Anaconda’s cumulative artifact total, the
retained PyPI daily series excluding mirrors, observed BioContainers/Quay image
pulls in the reported stats window, cumulative Zenodo release-record downloads,
and available GitHub Release asset downloads. The retained PyPI series carries
forward dated observations from prior snapshots because the provider exposes a
finite history window; it begins with the first retained snapshot and is not a
provider lifetime total. It is deliberately not deduplicated across channels:
one install that retrieves from two channels counts as two provider events.
Because the Quay component is also a bounded observation window, it must not be
described as cumulative lifetime pulls.

GHCR remains tracked as a public distribution channel but is not added to the
aggregate because its unauthenticated endpoint does not expose a usable pull
count. BioContainers/Quay is counted when the public `includeStats=true`
endpoint is available. That value is the sum of daily image-pull counts in the
provider’s observed stats window, not a lifetime estimate; the snapshot keeps
the dates, public tag inventory, and explicit empty version/platform
breakdowns rather than guessing attribution.

Zenodo’s counter is recorded at the release-record level. The current record
contains one versioned source archive, so the snapshot can show its release
version and archive platform; the provider counter still does not establish
distinct users or runtime use.

GitHub stars and forks are retained separately as public discovery signals and
are never converted into package downloads or user estimates. If a provider is
temporarily unavailable, the snapshot records `availability: unavailable` and
the error rather than silently treating the missing response as zero.

The PyPI recent-month endpoint is also non-critical. If PyPI Stats rate-limits
that supplemental request, the snapshot preserves the retained no-mirror daily
series and records `recent_availability: unavailable` with the provider error;
it does not turn the missing recent value into zero.

## Evidence ladder

Use the download snapshot for distribution progress. Add separate, linkable
evidence before making stronger claims:

- an external issue, pull request, or workflow reference for discovery;
- a named, permissioned public use record for a real external run;
- repeat use, a citation, or an accepted workflow integration for adoption.

Do not add silent telemetry to the package to fill this gap. If a user chooses
to report a run, record the version, channel, workflow, platform, and outcome
without collecting sequencing data.

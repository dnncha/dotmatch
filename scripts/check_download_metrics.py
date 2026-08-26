#!/usr/bin/env python3
"""Collect public DotMatch distribution and discovery metrics.

The provider counts in this report are event counts, not distinct users. The
report deliberately keeps package downloads separate from weaker public
discovery signals such as GitHub stars.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import time
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable
from urllib.error import HTTPError, URLError


ANACONDA_PACKAGE_URL = "https://api.anaconda.org/package/bioconda/dotmatch"
ANACONDA_REPOCORE_URL = "https://api.anaconda.org/repocore/channels/bioconda/artifacts/conda/dotmatch"
PYPI_JSON_URL = "https://pypi.org/pypi/dotmatch/json"
PYPISTATS_RECENT_URL = "https://pypistats.org/api/packages/dotmatch/recent?period=month"
PYPISTATS_OVERALL_URL = "https://pypistats.org/api/packages/dotmatch/overall?mirrors=false"
BIOCONTAINERS_STATS_URL = "https://quay.io/api/v1/repository/biocontainers/dotmatch?includeStats=true"
ZENODO_RECORD_URL = "https://zenodo.org/api/records/21511337"
GITHUB_REPO_URL = "https://api.github.com/repos/dnncha/dotmatch"
GITHUB_RELEASES_URL = "https://api.github.com/repos/dnncha/dotmatch/releases?per_page=100"
FETCH_ATTEMPTS = 3
RETRYABLE_HTTP_STATUS = {408, 425, 429, 500, 502, 503, 504}

JsonFetcher = Callable[[str], Any]


def fetch_json(url: str) -> Any:
    headers = {
        "Accept": "application/vnd.github+json, application/json",
        "User-Agent": "DotMatch download metrics",
    }
    github_token = os.environ.get("GITHUB_TOKEN")
    if github_token and url.startswith("https://api.github.com/"):
        headers["Authorization"] = f"Bearer {github_token}"
    request = urllib.request.Request(
        url,
        headers=headers,
    )
    last_error: Exception | None = None
    for attempt in range(FETCH_ATTEMPTS):
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                return json.loads(response.read().decode("utf-8"))
        except HTTPError as exc:
            last_error = exc
            if exc.code not in RETRYABLE_HTTP_STATUS:
                raise
        except (URLError, TimeoutError) as exc:
            last_error = exc
        if attempt < FETCH_ATTEMPTS - 1:
            time.sleep(2**attempt)
    assert last_error is not None
    raise last_error


def fetch_optional_json(url: str) -> Any:
    """Fetch a non-critical discovery source without turning outage into zero."""

    try:
        return fetch_json(url)
    except Exception as exc:  # pragma: no cover - the exact provider error varies
        return {"_fetch_error": f"{type(exc).__name__}: {exc}"}


def integer(value: Any) -> int:
    try:
        return int(value or 0)
    except (TypeError, ValueError):
        return 0


def artifact_platform(basename: str) -> str:
    return basename.split("/", 1)[0] if "/" in basename else "unknown"


def artifact_python(basename: str) -> str:
    match = re.search(r"-py(\d{3})", basename)
    if match is None:
        return "not_applicable"
    digits = match.group(1)
    return f"{digits[0]}.{digits[1:]}"


def aggregate_files(files: list[dict[str, Any]]) -> dict[str, Any]:
    by_version: defaultdict[str, dict[str, int]] = defaultdict(lambda: {"downloads": 0, "artifact_files": 0})
    by_platform: defaultdict[str, dict[str, int]] = defaultdict(lambda: {"downloads": 0, "artifact_files": 0})
    by_python: defaultdict[str, dict[str, int]] = defaultdict(lambda: {"downloads": 0, "artifact_files": 0})
    artifacts: list[dict[str, Any]] = []

    for item in files:
        basename = str(item.get("basename") or "")
        version = str(item.get("version") or "unknown")
        platform = artifact_platform(basename)
        python_version = artifact_python(basename)
        downloads = integer(item.get("ndownloads"))
        by_version[version]["downloads"] += downloads
        by_version[version]["artifact_files"] += 1
        by_platform[platform]["downloads"] += downloads
        by_platform[platform]["artifact_files"] += 1
        by_python[python_version]["downloads"] += downloads
        by_python[python_version]["artifact_files"] += 1
        artifacts.append(
            {
                "version": version,
                "platform": platform,
                "python": python_version,
                "basename": basename,
                "downloads": downloads,
                "upload_time": item.get("upload_time"),
            }
        )

    return {
        "total_downloads": sum(item["downloads"] for item in artifacts),
        "artifact_files": len(artifacts),
        "by_version": dict(sorted(by_version.items(), reverse=True)),
        "by_platform": dict(sorted(by_platform.items())),
        "by_python": dict(sorted(by_python.items())),
        "artifacts": sorted(artifacts, key=lambda item: (item["version"], item["platform"], item["basename"])),
    }


def aggregate_pypi_downloads(payload: dict[str, Any]) -> int:
    """Sum PyPI Stats daily downloads while excluding mirror traffic."""

    return sum(item["downloads"] for item in pypi_daily_downloads(payload))


def pypi_daily_downloads(payload: dict[str, Any]) -> list[dict[str, Any]]:
    """Return the dated PyPI no-mirror series without dropping its time axis."""

    entries = payload.get("data") or []
    return [
        {
            "date": str(item.get("date") or "unknown"),
            "downloads": integer(item.get("downloads")),
        }
        for item in entries
        if isinstance(item, dict) and item.get("category") == "without_mirrors"
    ]


def merge_pypi_daily_downloads(
    current: list[dict[str, Any]], previous_snapshot: dict[str, Any] | None = None
) -> list[dict[str, Any]]:
    """Retain dated observations that have rolled out of PyPI Stats' window.

    PyPI Stats exposes a finite recent-history window rather than a lifetime
    counter.  The scheduled metrics job therefore carries forward dated
    observations from the previous public snapshot and lets the newest
    provider response replace overlapping dates.  This is an observed
    cumulative series beginning with the first retained snapshot, not a claim
    about PyPI's lifetime total.
    """

    previous_data = previous_snapshot if isinstance(previous_snapshot, dict) else {}
    previous_channel = ((previous_data.get("channels") or {}).get("pypi") or {})
    previous = previous_channel.get("daily_downloads_excluding_mirrors") or []
    dated: dict[str, dict[str, Any]] = {}
    unknown: list[dict[str, Any]] = []

    for item in previous:
        if not isinstance(item, dict):
            continue
        date = str(item.get("date") or "unknown")
        normalized = {"date": date, "downloads": integer(item.get("downloads"))}
        if date != "unknown":
            dated[date] = normalized

    for item in current:
        date = str(item.get("date") or "unknown")
        normalized = {"date": date, "downloads": integer(item.get("downloads"))}
        if date == "unknown":
            unknown.append(normalized)
        else:
            dated[date] = normalized

    return [*sorted(dated.values(), key=lambda item: item["date"]), *unknown]


def release_asset_platform(name: str) -> str:
    lowered = name.lower()
    if "linux" in lowered or "manylinux" in lowered or "musllinux" in lowered:
        return "linux"
    if "macos" in lowered or "darwin" in lowered:
        return "macos"
    if "win" in lowered:
        return "windows"
    if name.endswith(".whl"):
        return "wheel_other"
    if name.endswith(".tar.gz") or name.endswith(".zip"):
        return "source_or_archive"
    return "metadata_or_other"


def aggregate_github_release_assets(payload: Any) -> dict[str, Any]:
    """Aggregate public GitHub Release asset retrievals without counting drafts."""

    if not isinstance(payload, list):
        return {
            "availability": "unavailable",
            "fetch_error": payload.get("_fetch_error") if isinstance(payload, dict) else "invalid response",
            "total_downloads": None,
            "release_count": None,
            "asset_files": None,
            "by_version": {},
            "by_platform": {},
            "assets": [],
        }

    by_version: defaultdict[str, dict[str, int]] = defaultdict(lambda: {"downloads": 0, "asset_files": 0})
    by_platform: defaultdict[str, dict[str, int]] = defaultdict(lambda: {"downloads": 0, "asset_files": 0})
    assets: list[dict[str, Any]] = []
    release_count = 0

    for release in payload:
        if not isinstance(release, dict) or release.get("draft") or not release.get("published_at"):
            continue
        release_count += 1
        tag = str(release.get("tag_name") or "unknown")
        version = tag.removeprefix("v")
        for item in release.get("assets") or []:
            if not isinstance(item, dict):
                continue
            name = str(item.get("name") or "unknown")
            platform = release_asset_platform(name)
            downloads = integer(item.get("download_count"))
            by_version[version]["downloads"] += downloads
            by_version[version]["asset_files"] += 1
            by_platform[platform]["downloads"] += downloads
            by_platform[platform]["asset_files"] += 1
            assets.append(
                {
                    "version": version,
                    "platform": platform,
                    "name": name,
                    "downloads": downloads,
                    "release_url": release.get("html_url"),
                    "published_at": release.get("published_at"),
                }
            )

    return {
        "availability": "available",
        "fetch_error": None,
        "total_downloads": sum(item["downloads"] for item in assets),
        "release_count": release_count,
        "asset_files": len(assets),
        "by_version": dict(sorted(by_version.items(), reverse=True)),
        "by_platform": dict(sorted(by_platform.items())),
        "assets": sorted(assets, key=lambda item: (item["version"], item["name"])),
    }


def aggregate_biocontainers_stats(payload: Any) -> dict[str, Any]:
    """Aggregate Quay's public pull-count series without inventing tag attribution."""

    if not isinstance(payload, dict):
        return {
            "availability": "unavailable",
            "fetch_error": "invalid response",
            "total_downloads": None,
            "daily_downloads": [],
            "observed_start": None,
            "observed_end": None,
            "tag_count": None,
            "tag_names": [],
            "by_version": {},
            "by_platform": {},
        }

    provider_error = payload.get("_fetch_error")
    stats = payload.get("stats")
    if provider_error:
        return {
            "availability": "unavailable",
            "fetch_error": provider_error,
            "total_downloads": None,
            "daily_downloads": [],
            "observed_start": None,
            "observed_end": None,
            "tag_count": None,
            "tag_names": [],
            "by_version": {},
            "by_platform": {},
        }
    if not isinstance(stats, list):
        return {
            "availability": "unavailable",
            "fetch_error": "stats missing from provider response",
            "total_downloads": None,
            "daily_downloads": [],
            "observed_start": None,
            "observed_end": None,
            "tag_count": None,
            "tag_names": [],
            "by_version": {},
            "by_platform": {},
        }

    daily_downloads = [
        {
            "date": str(item.get("date") or "unknown"),
            "downloads": integer(item.get("count")),
        }
        for item in stats
        if isinstance(item, dict)
    ]
    observed_dates = sorted(
        item["date"] for item in daily_downloads if item["date"] != "unknown"
    )
    tags = payload.get("tags")
    tag_names = sorted(str(name) for name in tags) if isinstance(tags, dict) else []

    return {
        "availability": "available",
        "fetch_error": None,
        "total_downloads": sum(item["downloads"] for item in daily_downloads),
        "daily_downloads": daily_downloads,
        "observed_start": observed_dates[0] if observed_dates else None,
        "observed_end": observed_dates[-1] if observed_dates else None,
        "tag_count": len(tag_names),
        "tag_names": tag_names,
        "by_version": {},
        "by_platform": {},
    }


def aggregate_zenodo_stats(payload: Any) -> dict[str, Any]:
    """Aggregate the release-record download counter exposed by Zenodo."""

    if not isinstance(payload, dict):
        return {
            "availability": "unavailable",
            "fetch_error": "invalid response",
            "total_downloads": None,
            "record_id": None,
            "doi": None,
            "version": None,
            "file_count": None,
            "files": [],
            "by_version": {},
            "by_platform": {},
        }

    provider_error = payload.get("_fetch_error")
    stats = payload.get("stats")
    if provider_error:
        return {
            "availability": "unavailable",
            "fetch_error": provider_error,
            "total_downloads": None,
            "record_id": None,
            "doi": None,
            "version": None,
            "file_count": None,
            "files": [],
            "by_version": {},
            "by_platform": {},
        }
    if not isinstance(stats, dict):
        return {
            "availability": "unavailable",
            "fetch_error": "stats missing from provider response",
            "total_downloads": None,
            "record_id": None,
            "doi": None,
            "version": None,
            "file_count": None,
            "files": [],
            "by_version": {},
            "by_platform": {},
        }

    metadata = payload.get("metadata") or {}
    version = str(metadata.get("version") or "unknown")
    raw_files = payload.get("files") or []
    files = [
        {
            "name": str(item.get("key") or "unknown"),
            "size": integer(item.get("size")),
            "checksum": item.get("checksum"),
            "platform": release_asset_platform(str(item.get("key") or "")),
        }
        for item in raw_files
        if isinstance(item, dict)
    ]
    total_downloads = integer(stats.get("downloads"))
    by_platform = {}
    if len(files) == 1:
        platform = files[0]["platform"]
        by_platform[platform] = {"downloads": total_downloads, "artifact_files": 1}

    return {
        "availability": "available",
        "fetch_error": None,
        "total_downloads": total_downloads,
        "record_id": payload.get("id"),
        "doi": payload.get("doi"),
        "version": version,
        "file_count": len(files),
        "files": files,
        "provider_unique_downloads": integer(stats.get("unique_downloads")),
        "provider_views": integer(stats.get("views")),
        "by_version": {
            version: {"downloads": total_downloads, "artifact_files": len(files)}
        },
        "by_platform": by_platform,
    }


def collect_snapshot(
    fetch: JsonFetcher | None = None,
    previous_snapshot: dict[str, Any] | None = None,
) -> dict[str, Any]:
    fetch = fetch_json if fetch is None else fetch
    anaconda_files = fetch(ANACONDA_PACKAGE_URL).get("files") or []
    anaconda_package = fetch(ANACONDA_REPOCORE_URL)
    pypi = fetch(PYPI_JSON_URL)
    pypistats_recent = fetch_optional_json(PYPISTATS_RECENT_URL)
    pypistats_overall = fetch(PYPISTATS_OVERALL_URL)
    biocontainers = fetch_optional_json(BIOCONTAINERS_STATS_URL)
    zenodo = fetch_optional_json(ZENODO_RECORD_URL)
    github = fetch_optional_json(GITHUB_REPO_URL)
    github_releases = fetch_optional_json(GITHUB_RELEASES_URL)
    anaconda_breakdown = aggregate_files(anaconda_files)
    biocontainers_breakdown = aggregate_biocontainers_stats(biocontainers)
    zenodo_breakdown = aggregate_zenodo_stats(zenodo)
    github_release_breakdown = aggregate_github_release_assets(github_releases)
    anaconda_downloads = integer(anaconda_package.get("download_count"))
    pypi_current_daily = pypi_daily_downloads(pypistats_overall)
    pypi_daily = merge_pypi_daily_downloads(pypi_current_daily, previous_snapshot)
    pypi_current_window_downloads = sum(item["downloads"] for item in pypi_current_daily)
    pypi_downloads = sum(item["downloads"] for item in pypi_daily)
    pypi_recent_error = pypistats_recent.get("_fetch_error") if isinstance(pypistats_recent, dict) else "invalid response"
    github_release_downloads = github_release_breakdown["total_downloads"] or 0
    biocontainers_downloads = biocontainers_breakdown["total_downloads"] or 0
    zenodo_downloads = zenodo_breakdown["total_downloads"] or 0
    aggregate_downloads = (
        anaconda_downloads
        + pypi_downloads
        + biocontainers_downloads
        + zenodo_downloads
        + github_release_downloads
    )
    counted_channels = ["anaconda_bioconda", "pypi"]
    if biocontainers_breakdown["availability"] == "available":
        counted_channels.append("biocontainers")
    if zenodo_breakdown["availability"] == "available":
        counted_channels.append("zenodo")
    if github_release_breakdown["availability"] == "available":
        counted_channels.append("github_releases")

    return {
        "schema_version": 2,
        "captured_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "aggregate_downloads": {
            "scope": "aggregate public package download events; not distinct users",
            "reported_downloads": aggregate_downloads,
            "counted_channels": counted_channels,
            "aggregation_note": (
                "Sums provider totals without cross-channel deduplication; retrievals from "
                "multiple channels remain separate provider events. The BioContainers value "
                "is the observed Quay stats-window total, not a lifetime estimate. Zenodo "
                "is counted from its cumulative release-record counter."
            ),
        },
        "channels": {
            "anaconda_bioconda": {
                "metric": "cumulative artifact downloads",
                "classification": "mixed_unknown",
                "source_url": ANACONDA_PACKAGE_URL,
                "repocore_source_url": ANACONDA_REPOCORE_URL,
                "latest_version": anaconda_package.get("latest_version"),
                "metadata_updated_at": anaconda_package.get("updated_at"),
                "download_count_reported_by_provider": integer(anaconda_package.get("download_count")),
                "breakdown": anaconda_breakdown,
                "automation_indicators": [
                    "counts are split across package artifacts, platforms, and Python builds",
                    "the provider does not expose client identity or runtime use",
                ],
                "human_discovery_indicators": [],
            },
            "pypi": {
                "metric": "observed cumulative downloads excluding mirrors from retained daily series",
                "classification": "mixed_unknown",
                "source_url": PYPISTATS_OVERALL_URL,
                "recent_source_url": PYPISTATS_RECENT_URL,
                "latest_version": (pypi.get("info") or {}).get("version"),
                "cumulative_downloads_excluding_mirrors": pypi_downloads,
                "current_window_downloads_excluding_mirrors": pypi_current_window_downloads,
                "daily_downloads_excluding_mirrors": pypi_daily,
                "series_scope": (
                    "dated observations retained from the public metrics snapshot; starts at the "
                    "first retained snapshot and is not a provider lifetime counter"
                ),
                "recent_availability": "unavailable" if pypi_recent_error else "available",
                "recent_month_downloads": (
                    None if pypi_recent_error else integer((pypistats_recent.get("data") or {}).get("last_month"))
                ),
                "recent_fetch_error": pypi_recent_error,
                "automation_indicators": [
                    "the public recent endpoint does not identify clients or workflows",
                    "the retained series explicitly excludes mirrors but does not identify clients or workflows",
                ],
                "human_discovery_indicators": [],
            },
            "github_releases": {
                "metric": "public GitHub Release asset downloads",
                "classification": "mixed_unknown",
                "source_url": GITHUB_RELEASES_URL,
                "download_count_reported_by_provider": github_release_breakdown["total_downloads"],
                "included_in_aggregate_total": github_release_breakdown["availability"] == "available",
                "availability": github_release_breakdown["availability"],
                "fetch_error": github_release_breakdown["fetch_error"],
                "breakdown": github_release_breakdown,
                "automation_indicators": [
                    "release assets may be retrieved by CI, mirrors, or direct users",
                    "GitHub does not expose client identity or runtime use in this count",
                ],
                "human_discovery_indicators": [],
            },
            "ghcr": {
                "metric": "public image availability; pull count not exposed by the unauthenticated registry endpoint",
                "classification": "distribution_channel_unmeasured",
                "source_url": "https://github.com/dnncha/dotmatch/pkgs/container/dotmatch",
                "count_source_url": "https://api.github.com/users/dnncha/packages/container/dotmatch",
                "count_status": "authentication_required",
                "included_in_aggregate_total": False,
            },
            "biocontainers": {
                "metric": "observed container pull events reported by Quay's repository stats window",
                "classification": "mixed_unknown",
                "source_url": "https://quay.io/repository/biocontainers/dotmatch",
                "count_source_url": BIOCONTAINERS_STATS_URL,
                "download_count_reported_by_provider": biocontainers_breakdown["total_downloads"],
                "included_in_aggregate_total": biocontainers_breakdown["availability"] == "available",
                "availability": biocontainers_breakdown["availability"],
                "fetch_error": biocontainers_breakdown["fetch_error"],
                "breakdown": biocontainers_breakdown,
                "automation_indicators": [
                    "pull counts cover all public tags in the provider's reported stats window",
                    (
                        "the provider does not expose client identity, runtime use, version "
                        "attribution, or platform attribution"
                    ),
                ],
                "human_discovery_indicators": [],
            },
            "zenodo": {
                "metric": "cumulative release-record downloads",
                "classification": "mixed_unknown",
                "source_url": "https://doi.org/10.5281/zenodo.21511337",
                "count_source_url": ZENODO_RECORD_URL,
                "download_count_reported_by_provider": zenodo_breakdown["total_downloads"],
                "included_in_aggregate_total": zenodo_breakdown["availability"] == "available",
                "availability": zenodo_breakdown["availability"],
                "fetch_error": zenodo_breakdown["fetch_error"],
                "breakdown": zenodo_breakdown,
                "automation_indicators": [
                    "release archives may be retrieved by CI, mirrors, or direct users",
                    "Zenodo's counter does not establish client identity or runtime use",
                ],
                "human_discovery_indicators": [
                    "record views are retained as discovery context, not download evidence"
                ],
            },
        },
        "discovery_signals": {
            "github": {
                "source_url": GITHUB_REPO_URL,
                "availability": "unavailable" if github.get("_fetch_error") else "available",
                "stars": integer(github.get("stargazers_count")),
                "forks": integer(github.get("forks_count")),
                "classification": "weak_public_discovery_signal",
                "fetch_error": github.get("_fetch_error"),
            },
        },
        "interpretation": {
            "download_counts_are": "provider-reported package retrieval events",
            "aggregate_total_is": (
                "the sum of Anaconda cumulative artifact downloads, the retained PyPI daily series "
                "excluding mirrors, observed BioContainers image pulls in Quay's stats window, "
                "cumulative Zenodo release-record downloads, and available GitHub Release "
                "asset downloads"
            ),
            "download_counts_are_not": [
                "unique users",
                "unique organizations",
                "page views",
                "proof of scientific use",
            ],
            "human_use_requires_separate_evidence": True,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        help="write the JSON snapshot to this path instead of stdout",
    )
    args = parser.parse_args()
    previous_snapshot = None
    if args.output is not None and args.output.exists():
        try:
            previous_snapshot = json.loads(args.output.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous_snapshot = None
    snapshot = collect_snapshot(previous_snapshot=previous_snapshot)
    encoded = json.dumps(snapshot, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
        print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Check required platforms in an OCI image index or Docker manifest list."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def available_platforms(manifest: object) -> set[str]:
    if not isinstance(manifest, dict):
        raise ValueError("OCI manifest must be a JSON object")
    if int(manifest.get("schemaVersion") or 0) != 2:
        raise ValueError("OCI manifest must use schemaVersion 2")
    descriptors = manifest.get("manifests")
    if not isinstance(descriptors, list):
        raise ValueError("OCI manifest must be an image index or manifest list")

    platforms: set[str] = set()
    for descriptor in descriptors:
        if not isinstance(descriptor, dict):
            continue
        platform = descriptor.get("platform")
        if not isinstance(platform, dict):
            continue
        operating_system = str(platform.get("os") or "").strip()
        architecture = str(platform.get("architecture") or "").strip()
        if operating_system and architecture:
            platforms.add(f"{operating_system}/{architecture}")
    return platforms


def check_manifest(manifest: object, required_platforms: list[str]) -> set[str]:
    platforms = available_platforms(manifest)
    missing = sorted(set(required_platforms) - platforms)
    if missing:
        raise ValueError(f"OCI manifest is missing required platform(s): {', '.join(missing)}")
    return platforms


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="OCI image index or Docker manifest-list JSON file")
    parser.add_argument(
        "--require-platform",
        action="append",
        dest="required_platforms",
        default=[],
        metavar="OS/ARCH",
        help="required platform, for example linux/arm64; may be repeated",
    )
    args = parser.parse_args()

    try:
        manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
        platforms = check_manifest(manifest, args.required_platforms)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"OCI MANIFEST: FAIL ({exc})")
        return 1

    rendered = ", ".join(sorted(platforms)) or "none"
    print(f"OCI MANIFEST: PASS ({rendered})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

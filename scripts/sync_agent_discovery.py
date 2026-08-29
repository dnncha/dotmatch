#!/usr/bin/env python3
"""Copy the canonical agent discovery files to each published surface."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COPY_TARGETS = {
    "llms.txt": ["docs/llms.txt", "public/llms.txt"],
    "llms-full.txt": ["docs/llms-full.txt", "public/llms-full.txt"],
    "agent-capabilities.json": [
        "docs/agent-capabilities.json",
        "public/agent-capabilities.json",
        "python/dotmatch/data/agent-capabilities.json",
    ],
    "agent-capabilities.schema.json": [
        "docs/agent-capabilities.schema.json",
        "public/agent-capabilities.schema.json",
        "python/dotmatch/data/agent-capabilities.schema.json",
    ],
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail instead of updating stale copies")
    args = parser.parse_args()

    stale: list[str] = []
    for source_name, target_names in COPY_TARGETS.items():
        source = ROOT / source_name
        source_bytes = source.read_bytes()
        for target_name in target_names:
            target = ROOT / target_name
            if target.is_file() and target.read_bytes() == source_bytes:
                continue
            stale.append(target_name)
            if not args.check:
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, target)

    if args.check and stale:
        print("Agent discovery copies are stale: " + ", ".join(stale))
        print("Run: python3 scripts/sync_agent_discovery.py")
        return 1

    if stale:
        print("Updated agent discovery copies: " + ", ".join(stale))
    else:
        print("Agent discovery copies are aligned")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

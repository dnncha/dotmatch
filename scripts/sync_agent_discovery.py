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
    "agent-capabilities.v1.json": [
        "docs/agent-capabilities.v1.json",
        "public/agent-capabilities.v1.json",
    ],
    "agent-capabilities.v1.schema.json": [
        "docs/agent-capabilities.v1.schema.json",
        "public/agent-capabilities.v1.schema.json",
    ],
    "agent-tools.json": [
        "docs/agent-tools.json",
        "public/agent-tools.json",
        "python/dotmatch/data/agent-tools.json",
    ],
    "agent-tools.schema.json": [
        "docs/agent-tools.schema.json",
        "public/agent-tools.schema.json",
        "python/dotmatch/data/agent-tools.schema.json",
    ],
    "agent-reference-crispr.json": [
        "docs/agent-reference-crispr.json",
        "public/agent-reference-crispr.json",
    ],
}
SKILL_SOURCE = "extensions/codex/dotmatch-agent"
SKILL_TARGET = "python/dotmatch/data/codex-skill"
CLAUDE_SKILL_TARGETS = [
    "extensions/claude-code/dotmatch-agent",
    "python/dotmatch/data/claude-code-skill",
]


def _skill_files(root: Path, relative: str) -> dict[str, bytes]:
    directory = root / relative
    if not directory.is_dir():
        return {}
    return {
        path.relative_to(directory).as_posix(): path.read_bytes()
        for path in directory.rglob("*")
        if path.is_file()
    }


def _claude_skill_files(root: Path) -> dict[str, bytes]:
    return {
        name: content
        for name, content in _skill_files(root, SKILL_SOURCE).items()
        if name != "agents/openai.yaml"
    }


def _write_files(root: Path, relative: str, content: dict[str, bytes]) -> None:
    target = root / relative
    if target.exists():
        shutil.rmtree(target)
    for name, value in content.items():
        destination = target / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(value)


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

    expected_skill = _skill_files(ROOT, SKILL_SOURCE)
    installed_skill = _skill_files(ROOT, SKILL_TARGET)
    if expected_skill != installed_skill:
        stale.append(SKILL_TARGET)
        if not args.check:
            target = ROOT / SKILL_TARGET
            if target.exists():
                shutil.rmtree(target)
            shutil.copytree(ROOT / SKILL_SOURCE, target)

    expected_claude_skill = _claude_skill_files(ROOT)
    for target_name in CLAUDE_SKILL_TARGETS:
        if _skill_files(ROOT, target_name) == expected_claude_skill:
            continue
        stale.append(target_name)
        if not args.check:
            _write_files(ROOT, target_name, expected_claude_skill)

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

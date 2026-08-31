#!/usr/bin/env python3
"""Validate the DotMatch Codex skill and its packaged export copy."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "extensions" / "codex" / "dotmatch-agent"
PACKAGED = ROOT / "python" / "dotmatch" / "data" / "codex-skill"
CLAUDE_SKILL = ROOT / "extensions" / "claude-code" / "dotmatch-agent"
PACKAGED_CLAUDE = ROOT / "python" / "dotmatch" / "data" / "claude-code-skill"


def files(root: Path) -> dict[str, bytes]:
    return {path.relative_to(root).as_posix(): path.read_bytes() for path in root.rglob("*") if path.is_file()}


def main() -> int:
    failures: list[str] = []
    required = {
        "SKILL.md",
        "agents/openai.yaml",
        "references/crispr.md",
        "references/perturb-seq.md",
        "references/evidence-policy.md",
    }
    observed = set(files(SKILL)) if SKILL.is_dir() else set()
    if observed != required:
        failures.append(f"skill files must be exactly {sorted(required)}; observed {sorted(observed)}")
    if SKILL.is_dir() and PACKAGED.is_dir() and files(SKILL) != files(PACKAGED):
        failures.append("installed-package skill copy differs from the extension")
    elif not PACKAGED.is_dir():
        failures.append("installed-package skill copy is missing")

    claude_required = required - {"agents/openai.yaml"}
    claude_observed = set(files(CLAUDE_SKILL)) if CLAUDE_SKILL.is_dir() else set()
    if claude_observed != claude_required:
        failures.append(f"Claude Code skill files must be exactly {sorted(claude_required)}; observed {sorted(claude_observed)}")
    if CLAUDE_SKILL.is_dir() and PACKAGED_CLAUDE.is_dir() and files(CLAUDE_SKILL) != files(PACKAGED_CLAUDE):
        failures.append("installed-package Claude Code skill copy differs from the extension")
    elif not PACKAGED_CLAUDE.is_dir():
        failures.append("installed-package Claude Code skill copy is missing")

    skill_text = (SKILL / "SKILL.md").read_text(encoding="utf-8") if (SKILL / "SKILL.md").is_file() else ""
    for phrase in [
        "name: dotmatch-agent",
        "dotmatch agent tools --json",
        "Never edit target sequences",
        "references/crispr.md",
        "references/perturb-seq.md",
        "references/evidence-policy.md",
    ]:
        if phrase not in skill_text:
            failures.append(f"SKILL.md missing: {phrase}")
    yaml_text = (SKILL / "agents" / "openai.yaml").read_text(encoding="utf-8") if (SKILL / "agents" / "openai.yaml").is_file() else ""
    for phrase in ["display_name:", "short_description:", "default_prompt:", "$dotmatch-agent", "allow_implicit_invocation: true"]:
        if phrase not in yaml_text:
            failures.append(f"agents/openai.yaml missing: {phrase}")

    codex_root = Path(os.environ.get("CODEX_HOME", str(Path.home() / ".codex")))
    standard_validator = codex_root / "skills" / ".system" / "skill-creator" / "scripts" / "quick_validate.py"
    if standard_validator.is_file():
        for label, skill_path in (("Codex", SKILL), ("Claude Code", CLAUDE_SKILL)):
            result = subprocess.run(
                [sys.executable, str(standard_validator), str(skill_path)],
                cwd=ROOT,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            if result.returncode:
                failures.append(f"standard skill validator failed for {label}: " + result.stdout.strip())
            else:
                print(f"{label}: {result.stdout.strip()}")
    else:
        print("Standard Codex skill validator is unavailable on this host; repository checks still ran")

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    print("Codex and Claude Code skill structures, metadata, references, and packaged copies passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

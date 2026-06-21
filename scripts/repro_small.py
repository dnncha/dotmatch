#!/usr/bin/env python3
"""Build a compact reviewer reproducibility packet from local evidence files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


DEFAULT_COMMANDS = [
    "make test",
    "make cli-test",
    "make scientific-readiness-ready",
    "make evidence-gallery-ready",
    "make assay-evidence-ready",
]

EVIDENCE_PATHS = [
    "docs/assay-evidence.json",
    "docs/scientific-readiness.json",
    "docs/workflow-adoption.json",
    "docs/benchmarks/native/README.md",
    "docs/benchmarks/public_crispr/README.md",
    "docs/benchmarks/crispr_comparison/README.md",
    "docs/benchmarks/barcode_demux/README.md",
    "docs/evidence-gallery/README.md",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def run_command(command: str, root: Path) -> dict[str, object]:
    started = datetime.now(timezone.utc)
    proc = subprocess.run(
        command,
        cwd=root,
        shell=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    return {
        "command": command,
        "exit_code": proc.returncode,
        "started_utc": started.isoformat(timespec="seconds"),
        "output_tail": proc.stdout.splitlines()[-80:],
    }


def assay_rows(root: Path) -> list[dict[str, str]]:
    manifest = read_json(root / "docs" / "assay-evidence.json")
    rows: list[dict[str, str]] = []
    for assay in manifest.get("assays", []):
        if not isinstance(assay, dict):
            continue
        rows.append(
            {
                "assay": str(assay.get("name") or assay.get("id") or ""),
                "status": str(assay.get("status") or ""),
                "gate": ", ".join(str(item) for item in assay.get("gates", []) or []),
                "claim_boundary": str(assay.get("claim_boundary") or ""),
                "next_public_evidence": str(assay.get("next_public_evidence") or ""),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fields = ["assay", "status", "gate", "claim_boundary", "next_public_evidence"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(packet: dict[str, object], rows: list[dict[str, str]], path: Path) -> None:
    lines = [
        "# DotMatch Reviewer Reproducibility Packet",
        "",
        f"Generated UTC: `{packet['generated_utc']}`",
        f"Git commit: `{packet['git_commit']}`",
        f"Dirty worktree observed: `{packet['git_dirty']}`",
        "",
        "## Focused Reproduction",
        "",
        "Run the compact reviewer target:",
        "",
        "```bash",
        "make repro-small",
        "```",
        "",
        "This target builds the native CLI, runs native and CLI fixture tests,",
        "checks scientific-readiness guardrails, and validates the",
        "evidence-gallery and assay-evidence manifests. Public-data comparison",
        "gates remain separate because they depend on larger cached benchmark",
        "artifacts and external competitor runs.",
        "",
        "## CI Artifact",
        "",
        "The `ci` workflow uploads the packet as the `reviewer-repro-packet` artifact",
        "from the `reviewer-repro` job. The artifact is intentionally small enough for",
        "reviewers to inspect without rerunning full public-data benchmarks.",
        "",
        "## Assay And Resubmission Matrix",
        "",
        "| Assay | Status | Gate | Next public evidence |",
        "| --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            "| {assay} | {status} | {gate} | {next_public_evidence} |".format(
                assay=row["assay"].replace("|", "\\|"),
                status=row["status"].replace("|", "\\|"),
                gate=row["gate"].replace("|", "\\|"),
                next_public_evidence=row["next_public_evidence"].replace("|", "\\|"),
            )
        )
    lines.extend(
        [
            "",
            "## Evidence File Checksums",
            "",
            "| Path | SHA-256 |",
            "| --- | --- |",
        ]
    )
    for item in packet["evidence_files"]:
        lines.append(f"| `{item['path']}` | `{item['sha256']}` |")
    lines.extend(
        [
            "",
            "## Commands",
            "",
            "| Command | Exit code |",
            "| --- | --- |",
        ]
    )
    for item in packet["commands"]:
        exit_code = item.get("exit_code", "not-run")
        lines.append(f"| `{item['command']}` | `{exit_code}` |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="repository root")
    parser.add_argument("--out-dir", default="repro/small", help="packet output directory")
    parser.add_argument("--skip-run", action="store_true", help="write packet without running commands")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    out_dir = (root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    git_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    ).stdout.strip()
    git_dirty = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    ).stdout.strip()

    commands = []
    for command in DEFAULT_COMMANDS:
        if args.skip_run:
            commands.append({"command": command, "exit_code": "not-run"})
            continue
        result = run_command(command, root)
        commands.append(result)
        (out_dir / f"{command.replace(' ', '_').replace('/', '_')}.log").write_text(
            "\n".join(str(line) for line in result["output_tail"]) + "\n",
            encoding="utf-8",
        )
        if result["exit_code"] != 0:
            break

    evidence_files = []
    for relative in EVIDENCE_PATHS:
        path = root / relative
        if path.exists():
            evidence_files.append({"path": relative, "sha256": sha256(path), "bytes": path.stat().st_size})

    rows = assay_rows(root)
    write_csv(out_dir / "resubmission-matrix.csv", rows)

    packet: dict[str, object] = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": git_commit,
        "git_dirty": bool(git_dirty),
        "commands": commands,
        "evidence_files": evidence_files,
        "adoption_status": read_json(root / "docs" / "workflow-adoption.json"),
    }
    packet_text = json.dumps(packet, indent=2) + "\n"
    (out_dir / "reviewer-repro-packet.json").write_text(packet_text, encoding="utf-8")
    (out_dir / "repro_manifest.json").write_text(packet_text, encoding="utf-8")
    write_markdown(packet, rows, out_dir / "README.md")

    failing = [item for item in commands if isinstance(item.get("exit_code"), int) and item["exit_code"] != 0]
    if failing:
        print(f"Reviewer repro packet written to {out_dir}", file=sys.stderr)
        print(f"Failing command: {failing[0]['command']}", file=sys.stderr)
        return int(failing[0]["exit_code"])
    print(f"Reviewer repro packet written to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

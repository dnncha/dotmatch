#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def require(text: str, needle: str, label: str, failures: list[str]) -> None:
    if " ".join(needle.split()) not in " ".join(text.split()):
        failures.append(f"{label} missing: {needle}")


def main() -> int:
    failures: list[str] = []
    docs = read("docs/workbench.md")
    readme = read("README.md")
    package_json = read("apps/workbench/package.json")
    tauri_config = read("apps/workbench/src-tauri/tauri.conf.json")
    rust = read("apps/workbench/src-tauri/src/lib.rs")
    app = read("apps/workbench/src/App.tsx")

    for needle in [
        "optional local desktop app",
        "All sequencing data stays on the user's machine.",
        "The Workbench is not part of the Bioconda recipe",
        "DOTMATCH_WORKBENCH_DOTMATCH",
        "No hosted uploads, accounts, telemetry, cloud storage, or external workflow adoption claims are required",
        "workspace confinement",
    ]:
        require(docs, needle, "docs/workbench.md", failures)

    require(readme, "Optional local Workbench", "README.md", failures)
    require(readme, "separate from the Bioconda recipe", "README.md", failures)
    require(package_json, '"tauri:build"', "apps/workbench/package.json", failures)
    require(tauri_config, '"DotMatch Workbench"', "Tauri config", failures)
    require(rust, "validate_workspace_args", "Workbench backend", failures)
    require(rust, "build_dotmatch_command", "Workbench backend", failures)
    require(rust, "DOTMATCH_WORKBENCH_DOTMATCH", "Workbench backend", failures)
    require(app, "Data stays local", "Workbench UI", failures)
    require(app, "AssayInfer", "Workbench UI", failures)
    require(app, "AssayRun", "Workbench UI", failures)
    require(app, "AssayAutopsy", "Workbench UI", failures)

    combined_public = "\n".join([docs, readme])
    banned = [
        "available on Bioconda",
        "Bioconda package is available",
        "ToolShed",
        "nf-core accepted",
        "should not claim",
        "launch path",
    ]
    for phrase in banned:
        if phrase.lower() in combined_public.lower():
            failures.append(f"public Workbench copy contains banned phrase: {phrase}")
    if " ai " in f" {combined_public.lower()} ":
        failures.append("public Workbench copy contains standalone AI wording")

    if failures:
        print("Workbench surface check failed:")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print("Workbench surface check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

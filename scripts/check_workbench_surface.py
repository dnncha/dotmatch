#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(path: str, failures: list[str]) -> str:
    full = ROOT / path
    if not full.is_file():
        failures.append(f"missing required Workbench file: {path}")
        return ""
    return full.read_text(encoding="utf-8")


def main() -> int:
    failures: list[str] = []
    for path in [
        "docs/workbench.md",
        "apps/workbench/package.json",
        "apps/workbench/src-tauri/tauri.conf.json",
        "apps/workbench/src-tauri/Cargo.toml",
        "apps/workbench/src-tauri/src/lib.rs",
        "apps/workbench/src/App.tsx",
        "apps/workbench/src/lib/assayModel.ts",
        "apps/workbench/src/lib/results.ts",
        "apps/workbench/src/lib/workbenchApi.ts",
    ]:
        if not (ROOT / path).is_file():
            failures.append(f"missing required Workbench file: {path}")

    try:
        package_json = json.loads(read("apps/workbench/package.json", failures))
    except json.JSONDecodeError as exc:
        failures.append(f"apps/workbench/package.json is invalid JSON: {exc}")
        package_json = {}
    scripts = package_json.get("scripts", {})
    for script in ["build", "lint", "test", "tauri:build"]:
        if script not in scripts:
            failures.append(f"apps/workbench/package.json missing script: {script}")

    try:
        tauri_config = json.loads(read("apps/workbench/src-tauri/tauri.conf.json", failures))
    except json.JSONDecodeError as exc:
        failures.append(f"apps/workbench/src-tauri/tauri.conf.json is invalid JSON: {exc}")
        tauri_config = {}
    if tauri_config.get("productName") != "DotMatch Workbench":
        failures.append("Tauri config productName must be DotMatch Workbench")
    if not tauri_config.get("bundle", {}).get("targets"):
        failures.append("Tauri config must define bundle targets")

    rust = read("apps/workbench/src-tauri/src/lib.rs", failures)
    for symbol in [
        "canonical_workspace",
        "resolve_workspace_path",
        "validate_workspace_args",
        "build_dotmatch_command",
        "DOTMATCH_WORKBENCH_DOTMATCH",
        "run_workbench_command",
    ]:
        if symbol not in rust:
            failures.append(f"Workbench backend missing symbol: {symbol}")

    app = read("apps/workbench/src/App.tsx", failures)
    for symbol in ["AssayInfer", "AssayPlan", "AssayRun", "AssayAutopsy", "runWorkbenchCommand"]:
        if symbol not in app:
            failures.append(f"Workbench UI missing symbol: {symbol}")

    if failures:
        print("Workbench surface check failed:")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print("Workbench surface check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

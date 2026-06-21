#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    failures: list[str] = []

    if (ROOT / "apps" / "workbench").exists():
        failures.append("apps/workbench must live in the dotmatch-community repository")

    workbench_doc = ROOT / "docs" / "workbench.md"
    if not workbench_doc.is_file():
        failures.append("missing docs/workbench.md integration contract")
    else:
        text = workbench_doc.read_text(encoding="utf-8")
        for required in [
            "dotmatch-community",
            "DOTMATCH_WORKBENCH_DOTMATCH",
            "assay_manifest.json",
            "sample_qc.tsv",
        ]:
            if required not in text:
                failures.append(f"docs/workbench.md must mention {required}")

    readme = ROOT / "README.md"
    if readme.is_file():
        text = readme.read_text(encoding="utf-8")
        if "dotmatch-community" not in text:
            failures.append("README.md must point Workbench users to dotmatch-community")
    else:
        failures.append("missing README.md")

    if failures:
        print("Workbench boundary check failed:")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print("Workbench boundary check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

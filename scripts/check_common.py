"""Shared helpers for DotMatch repository gate scripts."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any


PLACEHOLDER_HOSTS = {"example.org", "example.com", "example.net"}


class AuditResult:
    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failures: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.failures


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_https_url(value: str) -> bool:
    return re.match(r"^https://[^ \t\r\n]+$", value) is not None


def uses_placeholder_host(value: str) -> bool:
    match = re.match(r"^https://([^/]+)", value)
    return bool(match and match.group(1).lower() in PLACEHOLDER_HOSTS)


def check_https_url(owner: str, field: str, value: str, result: AuditResult) -> bool:
    if not value or not is_https_url(value):
        result.failures.append(f"{owner} must declare {field} as an https URL")
        return False
    if uses_placeholder_host(value):
        result.failures.append(f"{owner} {field} must not use placeholder domains")
        return False
    return True


def make_targets(root: Path) -> set[str]:
    makefile = root / "Makefile"
    if not makefile.exists():
        return set()
    targets: set[str] = set()
    for line in makefile.read_text(encoding="utf-8").splitlines():
        match = re.match(r"^([A-Za-z0-9_.-]+):(?:\s|$)", line)
        if match:
            targets.add(match.group(1))
    return targets


def check_simple_make_target(root: Path, command: str, owner: str, result: AuditResult) -> None:
    parts = command.split()
    if len(parts) != 2 or parts[0] != "make":
        result.failures.append(f"{owner} must be a simple make target: {command}")
        return
    if parts[1] not in make_targets(root):
        result.failures.append(f"{owner} target is missing: {command}")

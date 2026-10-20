#!/usr/bin/env python3
"""Check runtime import boundaries without loading the package or native library."""

from __future__ import annotations

import ast
import sys
from importlib.util import resolve_name
from pathlib import Path
from typing import Iterator

ROOT = Path(__file__).resolve().parents[1]
# These modules must remain usable with only the standard library and the
# explicitly listed low-level dependencies, including inside function bodies.
LEAF_IMPORTS: dict[str, set[str]] = {
    "_bindings": set(),
    "_optional": set(),
    "_validation": set(),
    "fastq_io": set(),
    "target_io": set(),
    "count_io": {"_validation"},
    "native": set(),
}
HEAVY = {"anndata", "numpy", "pandas", "polars", "scipy", "multiqc"}
CLI = {"cli", "entrypoint", "__main__", "first_run"}
ORCHESTRATION = CLI | {"agent_tools"}
FORBIDDEN = {"quickdna", "assaycode", "scripts", "tests", "repos", "packages"}


def imports(
    tree: ast.AST, *, eager: bool = True
) -> Iterator[tuple[ast.Import | ast.ImportFrom, bool]]:
    """Function-body imports are lazy; imports inside module-level try/if are eager."""
    for node in ast.iter_child_nodes(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            yield node, eager
        else:
            lazy = isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda))
            yield from imports(node, eager=eager and not lazy)


def imported_names(node: ast.AST, module: str) -> list[str]:
    if isinstance(node, ast.Import):
        return [alias.name for alias in node.names]
    if isinstance(node, ast.ImportFrom):
        name = "." * node.level + (node.module or "")
        if node.level:
            name = resolve_name(name, module.rpartition(".")[0])
        if name == "dotmatch":
            return [f"dotmatch.{alias.name}" for alias in node.names]
        return [name]
    return []


def check(root: Path) -> list[str]:
    package = root / "python" / "dotmatch"
    paths = sorted(package.rglob("*.py"))
    if not paths:
        return ["python/dotmatch: no runtime modules found"]
    modules = {
        "dotmatch." + ".".join(path.relative_to(package).with_suffix("").parts): path
        for path in paths
    }
    graph: dict[str, set[str]] = {name: set() for name in modules}
    errors: list[str] = []
    for module, path in modules.items():
        short = module.removeprefix("dotmatch.")
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except SyntaxError as exc:
            errors.append(f"{path.relative_to(root)}:{exc.lineno}: {exc.msg}")
            continue
        for import_node, eager in imports(tree):
            for name in imported_names(import_node, module):
                top = name.split(".")[0]
                internal = name.removeprefix("dotmatch.") if top == "dotmatch" else None
                reason = ""
                if top in FORBIDDEN:
                    reason = (
                        "runtime must not depend on compatibility packages or repository tooling"
                    )
                elif (
                    internal is not None
                    and short in LEAF_IMPORTS
                    and internal not in LEAF_IMPORTS[short]
                ):
                    reason = "low-level module must not depend on higher-level runtime code"
                elif internal in CLI and short not in ORCHESTRATION:
                    reason = "library code must not depend on CLI orchestration"
                elif eager and top in HEAVY and short != "multiqc":
                    reason = "optional scientific integrations must be imported lazily"
                if reason:
                    errors.append(
                        f"{path.relative_to(root)}:{import_node.lineno}: {name}: {reason}"
                    )
                if eager and name in modules:
                    graph[module].add(name)
        # Loading native code has one owner. This also catches aliased imports.
        ctypes_aliases = {"ctypes"}
        loader_aliases: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                ctypes_aliases.update(
                    alias.asname or alias.name for alias in node.names if alias.name == "ctypes"
                )
            elif isinstance(node, ast.ImportFrom) and node.module == "ctypes":
                loader_aliases.update(
                    alias.asname or alias.name
                    for alias in node.names
                    if alias.name in {"CDLL", "PyDLL"}
                )
        if short != "_bindings":
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                fn = node.func
                direct = isinstance(fn, ast.Name) and fn.id in loader_aliases
                qualified = (
                    isinstance(fn, ast.Attribute)
                    and isinstance(fn.value, ast.Name)
                    and fn.value.id in ctypes_aliases
                    and fn.attr in {"CDLL", "PyDLL"}
                )
                if direct or qualified:
                    errors.append(
                        f"{path.relative_to(root)}:{node.lineno}: native library loading belongs in _bindings.py"
                    )

    # Detect eager cycles once, with a useful path, rather than executing imports.
    visited: set[str] = set()
    active: list[str] = []

    def visit(module: str) -> None:
        if module in active:
            cycle = active[active.index(module) :] + [module]
            errors.append("eager import cycle: " + " -> ".join(cycle))
            return
        if module in visited:
            return
        active.append(module)
        for dependency in sorted(graph[module]):
            visit(dependency)
        active.pop()
        visited.add(module)

    for module in sorted(graph):
        visit(module)
    return errors


def main() -> int:
    errors = check(ROOT)
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("Runtime architecture boundaries passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

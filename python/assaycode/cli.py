"""Compatibility-safe AssayCode command surface."""

from __future__ import annotations

import sys
from typing import Sequence

from dotmatch import __version__
from dotmatch import cli as _engine_cli


_SHORTCUTS = {"new", "infer", "check", "plan", "run", "start"}


def print_help() -> None:
    print(
        f"""AssayCode {__version__} — powered by the DotMatch engine

Compile, validate, decode, and diagnose sequencing assays built from known
guides, barcodes, primers, feature tags, and panel targets.

Usage:
  assaycode --help
  assaycode --version
  assaycode check assay.toml
  assaycode plan assay.toml
  assaycode run assay.toml
  assaycode start assay.toml
  assaycode assay <command> [options]
  assaycode engine <dotmatch-command> [options]

Assay workflow shortcuts:
  new       scaffold an AssayScript/AssaySpec project
  infer     infer a reviewable assay window from reads
  check     validate a specification and its referenced inputs
  plan      print the deterministic execution plan
  run       execute an already-reviewed specification
  start     check, plan, run, and write reliability artifacts

Specialized DotMatch workflows remain available unchanged:
  assaycode crispr ...
  assaycode barcode ...
  assaycode panel ...
  assaycode count ...
  assaycode demux ...

Compatibility:
  The dotmatch executable, Python package, native ABI, output schemas, DOI, and
  citation remain authoritative and supported. AssayCode is an additive
  assay-level identity, not a fork or destructive package rename.
"""
    )


def main(argv: Sequence[str] | None = None) -> int:
    raw_args = list(sys.argv[1:] if argv is None else argv)
    if not raw_args or raw_args[0] in {"-h", "--help", "help"}:
        print_help()
        return 0
    if raw_args == ["--version"]:
        print(f"assaycode {__version__} (DotMatch engine {__version__})")
        return 0
    if raw_args[0] == "engine":
        if len(raw_args) == 1:
            print("usage: assaycode engine <dotmatch-command> [options]", file=sys.stderr)
            return 2
        return _engine_cli.main(raw_args[1:])
    if raw_args[0] in _SHORTCUTS:
        return _engine_cli.main(["assay", *raw_args])
    return _engine_cli.main(raw_args)


if __name__ == "__main__":
    raise SystemExit(main())

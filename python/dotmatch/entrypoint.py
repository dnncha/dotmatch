"""Public CLI routing; leave the established analytical commands unchanged."""
from __future__ import annotations

import sys
from collections.abc import Sequence


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if args[:1] == ["sensitivity-review"]:
        from .sensitivity_review import main as review_main
        return review_main(args[1:])
    if args[:1] == ["demo"]:
        from .first_run import demo_main
        return demo_main(args[1:])
    if args[:1] == ["compare-counts"]:
        from .count_compare import main as compare_main
        return compare_main(args[1:])
    if args[:1] == ["evaluation-packet"]:
        from .evaluation_packet import main as evaluation_packet_main

        return evaluation_packet_main(args[1:])
    if args[:2] == ["crispr", "quickstart"]:
        from .first_run import quickstart_main
        return quickstart_main(args[2:])
    from .cli import main as existing_main
    result = existing_main(args)
    if args[:1] and args[0] in {"--help", "-h", "help"}:
        print("First run and workflow evaluation:\n"
              "  dotmatch demo --out-dir first-run/\n"
              "      Run a bundled synthetic example and check its expected results offline.\n"
              "  dotmatch compare-counts --help\n"
              "      Compare existing raw-count tables without replacing your workflow.\n"
              "  dotmatch evaluation-packet --help\n"
              "      Lock an independent comparison protocol before viewing outcomes.\n"
              "  dotmatch sensitivity-review --help\n"
              "      Rebuild an offline review from a completed sensitivity bundle.\n"
              "  dotmatch crispr quickstart --help\n"
              "      Prepare your own reads; --link-reads avoids copying large FASTQs.")
    return result

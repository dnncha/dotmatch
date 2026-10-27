"""Typed DotMatch quickstart.

Run from a checkout::

    DOTMATCH_LIB=$PWD/libdotmatch.dylib PYTHONPATH=$PWD/python python examples/python/typed_quickstart.py
"""

from __future__ import annotations

import json

import dotmatch as dm


def main() -> int:
    reads = ["ACGT", "ACGC", "TTTT"]
    targets = ["ACGT", "AGGT", "ACGA"]

    results = dm.assign(reads, targets, k=1)

    payload = {
        "version": dm.__version__,
        "statuses": [dm.status_name(result.status) for result in results],
        "typed": {
            "summary": {
                "total_reads": len(results),
                "unique": sum(1 for result in results if result.status == dm.MATCH_UNIQUE),
            }
        },
    }
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

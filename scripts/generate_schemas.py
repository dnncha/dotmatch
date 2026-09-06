"""Generate or check runtime JSON Schemas. Semantic validation is still authoritative."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from editwitness.cli import schema_for


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Fail instead of rewriting stale schemas")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1] / "src/editwitness/schemas"
    stale = []
    for kind in ("manifest", "analysis", "scan"):
        path = root / f"{kind}.schema.json"
        expected = schema_for(kind)
        if args.check:
            if not path.is_file() or json.loads(path.read_text(encoding="utf-8")) != expected:
                stale.append(str(path.relative_to(root.parent.parent.parent)))
        else:
            path.write_text(json.dumps(expected, indent=2) + "\n", encoding="utf-8")
    if stale:
        parser.exit(1, "Stale schemas: " + ", ".join(stale) + "\n")
    print("Schema snapshots match runtime contracts." if args.check else "Wrote schema snapshots.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

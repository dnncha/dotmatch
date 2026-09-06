"""Write checksums for fresh local build artifacts, without invoking any publisher."""
from __future__ import annotations

import hashlib
from pathlib import Path

root = Path(__file__).resolve().parents[1]
dist = root / "dist"
artifacts = sorted([*dist.glob("*.whl"), *dist.glob("*.tar.gz")])
if len(artifacts) != 2 or not any(p.suffix == ".whl" for p in artifacts):
    raise SystemExit("Expected exactly one wheel and one source distribution in a clean dist/.")
(dist / "SHA256SUMS").write_text("".join(
    f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n" for path in artifacts
), encoding="ascii")

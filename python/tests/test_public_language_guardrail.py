from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

PUBLIC_PATHS = [
    ROOT / "README.md",
    ROOT / "docs",
    ROOT / "app",
    ROOT / "public",
    ROOT / ".github",
]

PUBLIC_SUFFIXES = {".css", ".html", ".json", ".md", ".mdx", ".ts", ".tsx", ".yaml", ".yml"}
SKIP_PARTS = {"_build", "node_modules", "__pycache__"}

FORBIDDEN_PHRASES = [
    "adoption " + "evidence",
    "adoption " + "trust",
    "ai " + "slop",
    "best-in-class",
    "big " + "wins",
    "cutting-edge",
    "evidence-" + "bounded",
    "enterprise-grade",
    "future-proof",
    "game-changing",
    "industry " + "exposure",
    "industry-leading",
    "just works",
    "massive industry " + "impact",
    "next " + "wins",
    "pilot " + "conversations",
    "private " + "feedback",
    "quote-" + "approved",
    "revolutionary",
    "seamless",
    "turning private evaluation into public adoption " + "evidence",
    "without turning private feedback into public " + "evidence",
]


def _public_files():
    for path in PUBLIC_PATHS:
        if path.is_file():
            yield path
            continue
        if not path.exists():
            continue
        for candidate in path.rglob("*"):
            if not candidate.is_file() or candidate.suffix.lower() not in PUBLIC_SUFFIXES:
                continue
            if SKIP_PARTS.intersection(candidate.relative_to(ROOT).parts):
                continue
            yield candidate


def test_public_language_avoids_internal_process_phrasing():
    failures = []
    for path in _public_files():
        text = path.read_text(encoding="utf-8").lower()
        for phrase in FORBIDDEN_PHRASES:
            if phrase in text:
                failures.append(f"{path.relative_to(ROOT)} contains forbidden phrase: {phrase}")

    assert failures == []

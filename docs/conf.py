from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

project = "DotMatch"
author = "Donncha O'Toole"
copyright = "2026, Donncha O'Toole"
release = "0.3.0"
version = release

extensions = [
    "myst_parser",
]

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

master_doc = "index"

html_theme = "sphinx_rtd_theme"
html_title = "DotMatch documentation"
html_static_path = ["_static"]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]
myst_heading_anchors = 3

nitpicky = True
show_warning_types = True
suppress_warnings = [
    "myst.header",
    "myst.xref_missing",
]


NAVIGATION_DOCS = {
    "index",
    "getting-started",
    "command-reference",
    "tutorials/crispr-count-first-run",
    "tutorials/scverse-perturb-seq",
    "assayspec",
    "crispr-qc",
    "barcode-panel-design",
    "streaming-api",
    "schemas",
    "workbench",
    "trust-and-scope",
    "benchmarks/README",
    "methods-and-citation",
    "packaging",
}


def _mark_detail_page_as_orphan(app, docname, source):
    """Build detailed reports without placing them in the main user navigation."""
    if docname not in NAVIGATION_DOCS:
        source[0] = "---\norphan: true\n---\n\n" + source[0]


def setup(app):
    app.connect("source-read", _mark_detail_page_as_orphan)

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

project = "DotMatch"
author = "Donncha O'Toole"
copyright = "2026, Donncha O'Toole"
release = "0.1.8"
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
suppress_warnings = [
    "myst.header",
    "myst.xref_missing",
]

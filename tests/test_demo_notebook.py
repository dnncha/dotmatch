import json
import re
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_hosted_demo_installs_current_release() -> None:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    notebook = json.loads((ROOT / "demo.ipynb").read_text(encoding="utf-8"))
    source = "".join(
        line
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
        for line in cell["source"]
    )

    pins = re.findall(r'dotmatch==([0-9.]+)', source)
    assert pins == [project["project"]["version"]]

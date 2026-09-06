"""Run branch coverage including CLI subprocesses, retaining an auditable log."""
import os
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
env = dict(os.environ, PYTEST_DISABLE_PLUGIN_AUTOLOAD="1")
for args in (("erase",), ("run", "-m", "pytest", "-q"), ("combine",), ("report", "--fail-under=95"), ("json", "-o", "coverage.json")):
    subprocess.run([sys.executable, "-m", "coverage", *args], cwd=root, env=env, check=True)

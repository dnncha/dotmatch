import pathlib
import re
import sys


def requirement_list(pyproject_text: str, key: str) -> list[str]:
    block = re.search(rf"^{key} = \[(.*?)^\]", pyproject_text, re.S | re.M)
    if block is None:
        raise SystemExit(f"pyproject.toml has no {key} list")
    return [item.replace('\\"', '"') for item in re.findall(r'"((?:[^"\\]|\\.)+)"', block.group(1))]


def gate_requirements(pyproject_text: str) -> list[str]:
    return requirement_list(pyproject_text, "dependencies") + requirement_list(pyproject_text, "dev")


if __name__ == "__main__":
    root = pathlib.Path(__file__).resolve().parent.parent
    sys.stdout.write("\n".join(gate_requirements((root / "pyproject.toml").read_text(encoding="utf-8"))) + "\n")

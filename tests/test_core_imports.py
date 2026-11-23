"""The core imports neither the window toolkit nor the plotting library (spec invariant 6)."""

import pkgutil
import subprocess
import sys
from pathlib import Path

import pytest
from genomes import child_env

SRC = Path(__file__).resolve().parents[1] / "src"
NON_CORE = {"render", "cli", "gui"}
NOT_DEPENDENCIES = ("pandas", "sklearn", "scipy")
PROBE = (
    "import importlib, sys\n"
    "importlib.import_module(sys.argv[1])\n"
    "found = [name for name in sys.argv[2:] if name in sys.modules]\n"
    "print(','.join(found))\n"
    "sys.exit(1 if found else 0)\n"
)


def forbidden_for(module: str) -> tuple[str, ...]:
    if module not in NON_CORE:
        return ("tkinter", "matplotlib", *NOT_DEPENDENCIES)
    if module == "gui":
        return NOT_DEPENDENCIES
    return ("tkinter", *NOT_DEPENDENCIES)


def run_probe(module: str, env: dict[str, str]) -> subprocess.CompletedProcess:
    forbidden = forbidden_for(module.removeprefix("eulerase."))
    return subprocess.run(
        [sys.executable, "-c", PROBE, module, *forbidden],
        capture_output=True,
        text=True,
        timeout=60,
        env=env,
    )


MODULES = [m.name for m in pkgutil.iter_modules([str(SRC / "eulerase")])]


@pytest.mark.parametrize("module", ["eulerase", *(f"eulerase.{m}" for m in MODULES)])
def test_module_imports_nothing_it_must_not(module):
    result = run_probe(module, child_env())
    assert result.returncode == 0, f"{module}: {result.stdout.strip()} {result.stderr.strip()}"


def test_oracle_rejects_a_core_that_imports_matplotlib(tmp_path):
    package = tmp_path / "eulerase"
    package.mkdir()
    (package / "__init__.py").write_text("import matplotlib\n", encoding="utf-8")
    env = child_env()
    env["PYTHONPATH"] = str(tmp_path)
    env["MPLCONFIGDIR"] = str(tmp_path / "mpl")
    result = run_probe("eulerase", env)
    assert result.returncode == 1
    assert "matplotlib" in result.stdout.strip().split(",")

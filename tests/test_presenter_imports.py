"""The presenter imports no toolkit and no plotting library, and the probe can tell."""

import subprocess
import sys

import pytest
from genomes import child_env
from window_probe import blocked_code

SRC_PROBE = (
    "import importlib, sys\n"
    "importlib.import_module(sys.argv[1])\n"
    "found = [name for name in sys.argv[2:] if sys.modules.get(name) is not None]\n"
    "print(','.join(found))\n"
    "sys.exit(1 if found else 0)\n"
)
FORBIDDEN = ("matplotlib", "eulerase.render", "eulerase.gui.app", "pandas", "sklearn", "scipy")


def run_probe(module: str, env: dict[str, str]) -> subprocess.CompletedProcess:
    code = blocked_code(SRC_PROBE)
    return subprocess.run(
        [sys.executable, "-c", code, module, *FORBIDDEN],
        capture_output=True,
        text=True,
        timeout=60,
        env=env,
    )


def test_the_presenter_imports_no_toolkit_and_no_plotting_library():
    result = run_probe("eulerase.gui.presenter", child_env())
    assert result.returncode == 0, f"{result.stdout.strip()} {result.stderr.strip()}"


@pytest.mark.parametrize("library", ["matplotlib", "tkinter"])
def test_oracle_rejects_a_presenter_that_imports_a_toolkit_or_plotting_library(library, tmp_path):
    package = tmp_path / "eulerase"
    (package / "gui").mkdir(parents=True)
    (package / "__init__.py").write_text("", encoding="utf-8")
    (package / "gui" / "__init__.py").write_text("", encoding="utf-8")
    (package / "gui" / "presenter.py").write_text(f"import {library}\n", encoding="utf-8")
    env = child_env()
    env["PYTHONPATH"] = str(tmp_path)
    env["MPLCONFIGDIR"] = str(tmp_path / "mpl")
    result = run_probe("eulerase.gui.presenter", env)
    assert result.returncode != 0
    assert library in result.stderr + result.stdout

"""The window frame, checked without a window: import behaviour and source shape only."""

import ast
import inspect
import json
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("tkinter", reason="the window module imports tkinter; nothing here opens one")
import tkinter  # noqa: E402

assert getattr(tkinter.Tk, "eulerase_window_guard", False) and getattr(
    tkinter.Toplevel, "eulerase_window_guard", False
), "guard not installed"

from genomes import child_env  # noqa: E402
from window_probe import child_code  # noqa: E402

from eulerase import analysis  # noqa: E402
from eulerase.gui import app  # noqa: E402

APP_SOURCE = Path(app.__file__)
IMPORT_BODY = """
import importlib, json, threading
import matplotlib
backend = matplotlib.get_backend()
importlib.import_module("eulerase.gui.app")
print(json.dumps([threading.active_count(), matplotlib.get_backend() == backend]))
"""


def test_importing_the_window_creates_no_window_and_no_thread(tmp_path):
    env = child_env()
    env["MPLCONFIGDIR"] = str(tmp_path)
    env["MPLBACKEND"] = "Agg"
    done = subprocess.run(
        [sys.executable, "-c", child_code(IMPORT_BODY)],
        capture_output=True,
        text=True,
        timeout=120,
        env=env,
    )
    assert done.returncode != 3, "the window guard was not installed in the child"
    assert done.returncode == 0, done.stderr
    threads, same_backend = json.loads(done.stdout)
    assert threads == 1 and same_backend


def own_columns(source: str) -> list[str]:
    """String constants equal to a multi-word table header: tabs may only draw the frozen ones."""
    headers = {h for t in analysis.TABLES for h in analysis.columns(t) if " " in h}
    return sorted(
        {
            n.value
            for n in ast.walk(ast.parse(source))
            if isinstance(n, ast.Constant) and n.value in headers
        }
    )


def test_the_window_defines_no_table_columns_of_its_own():
    assert any(" " in h for t in analysis.TABLES for h in analysis.columns(t))
    assert own_columns(APP_SOURCE.read_text(encoding="utf-8")) == []


def test_oracle_rejects_a_module_with_its_own_columns():
    assert own_columns('COLUMNS = ("trail", "global step")') == ["global step"]


def test_the_window_gets_its_presenter_only_through_bind():
    assert list(inspect.signature(app.Window.__init__).parameters) == ["self", "root"]
    assert list(inspect.signature(app.Window.bind).parameters) == ["self", "presenter"]

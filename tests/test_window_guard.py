"""The window guard: no test, and no child interpreter, can open a real window."""

import subprocess
import sys

import pytest
from genomes import child_env
from window_probe import GUARD_PREAMBLE, MARKER_CHECK, child_code

tkinter = pytest.importorskip("tkinter", reason="the guard replaces tkinter classes")

NOT_INSTALLED = "the guard was not installed (child exit status 3)"


class RootWasBuilt(Exception):
    """Raised by the stand-in for the toolkit's root constructor: the guard let a call through."""


def refuse_to_build(*args, **kwargs):
    raise RootWasBuilt


def assert_guard_in_place():
    for name in ("Tk", "Toplevel"):
        assert getattr(getattr(tkinter, name), "eulerase_window_guard", False), NOT_INSTALLED


def run_child(code: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, timeout=60, env=child_env()
    )


def test_a_test_that_calls_tk_fails_at_once(monkeypatch):
    assert_guard_in_place()
    monkeypatch.setattr(tkinter._tkinter, "create", refuse_to_build)
    for name in ("Tk", "Toplevel"):
        with pytest.raises(AssertionError, match="a test tried to open a window"):
            getattr(tkinter, name)()


def test_the_dialog_modules_import_under_the_guard():
    """A fresh interpreter, so the imports run whatever the order of the tests."""
    body = "import tkinter.filedialog\nimport tkinter.simpledialog\n"
    result = run_child(child_code(body))
    assert result.returncode == 0, result.stderr


def test_a_test_can_replace_the_guard_with_its_own_fake(monkeypatch):
    class FakeRoot:
        pass

    assert_guard_in_place()
    monkeypatch.setattr(tkinter, "Tk", FakeRoot)
    assert isinstance(tkinter.Tk(), FakeRoot)
    monkeypatch.undo()
    assert_guard_in_place()


def test_the_marker_check_alone_stops_a_child_without_the_guard():
    result = run_child(MARKER_CHECK + "print('reached')\n")
    assert result.returncode == 3, result.stderr
    assert "reached" not in result.stdout


CALL_TK = """import tkinter

def refuse(*args, **kwargs):
    raise SystemExit(5)

tkinter._tkinter.create = refuse
try:
    tkinter.Tk()
except AssertionError as error:
    raise SystemExit(0 if "a test tried to open a window" in str(error) else 4)
raise SystemExit(4)
"""


def test_a_child_interpreter_guards_itself():
    result = run_child(GUARD_PREAMBLE + MARKER_CHECK + CALL_TK)
    assert result.returncode != 3, NOT_INSTALLED
    assert result.returncode != 5, "the guard did not raise: the root constructor was reached"
    assert result.returncode == 0, result.stderr

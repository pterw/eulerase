"""Suite-wide safety: no test opens a window, and matplotlib writes nothing under the home."""

import atexit
import os
import shutil
import tempfile

import pytest

os.environ["MPLBACKEND"] = "Agg"
_MPL_CONFIG_DIR = tempfile.mkdtemp(prefix="eulerase-mpl-")
os.environ["MPLCONFIGDIR"] = _MPL_CONFIG_DIR
atexit.register(shutil.rmtree, _MPL_CONFIG_DIR, ignore_errors=True)

GUARD_MESSAGE = "a test tried to open a window"


def _install_window_guard(patch=None) -> None:
    """Replace tkinter.Tk and tkinter.Toplevel with subclasses that raise before the parent runs."""
    try:
        import tkinter
    except ImportError:
        return

    def refuse(self, *args, **kwargs):
        raise AssertionError(GUARD_MESSAGE)

    for name in ("Tk", "Toplevel"):
        guard = type(
            name, (getattr(tkinter, name),), {"__init__": refuse, "eulerase_window_guard": True}
        )
        if patch is None:
            setattr(tkinter, name, guard)
        else:
            patch.setattr(tkinter, name, guard)


_install_window_guard()


@pytest.fixture(autouse=True)
def window_guard(monkeypatch):
    _install_window_guard(monkeypatch)

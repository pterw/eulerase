"""Test support: the code every child interpreter runs first, so none can open a window."""

GUARD_PREAMBLE = """\
try:
    import tkinter

    class _GuardedTk(tkinter.Tk):
        eulerase_window_guard = True

        def __init__(self, *args, **kwargs):
            raise AssertionError("a test tried to open a window")

    class _GuardedToplevel(tkinter.Toplevel):
        eulerase_window_guard = True

        def __init__(self, *args, **kwargs):
            raise AssertionError("a test tried to open a window")

    tkinter.Tk = _GuardedTk
    tkinter.Toplevel = _GuardedToplevel
except ImportError:
    pass
"""

MARKER_CHECK = """\
try:
    import tkinter
except ImportError:
    pass
else:
    import sys

    if not (
        getattr(tkinter.Tk, "eulerase_window_guard", False)
        and getattr(tkinter.Toplevel, "eulerase_window_guard", False)
    ):
        print("the window guard was not installed", file=sys.stderr)
        raise SystemExit(3)
"""


def child_code(body: str) -> str:
    """The guard, then the check that it is in place, then the body."""
    return f"{GUARD_PREAMBLE}{MARKER_CHECK}{body}"


def blocked_code(body: str) -> str:
    """For a probe that must see tkinter never imported: any import of it raises ImportError."""
    return f'import sys\nsys.modules["tkinter"] = None\n{body}'

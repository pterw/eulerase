"""Source rules for the window code, read from the syntax tree; nothing here imports the toolkit."""

import ast
import re
from collections.abc import Callable
from pathlib import Path

import pytest

GUI = Path(__file__).resolve().parents[1] / "src" / "eulerase" / "gui"
BARRED = re.compile(r"coverage|cov=|depth|bubble", re.IGNORECASE)
COLOUR = re.compile(r"#[0-9a-fA-F]{6}\b")
PRESENTER_BANNED = ("tkinter", "matplotlib", "eulerase.render", "eulerase.gui.app")
DYNAMIC_IMPORTS = ("__import__", "import_module")


def parse(source: str) -> ast.Module:
    return ast.parse(source)


def strings(source: str) -> list[str]:
    return [
        n.value
        for n in ast.walk(parse(source))
        if isinstance(n, ast.Constant) and isinstance(n.value, str)
    ]


def outside_the_basic_plane(source: str) -> list[str]:
    return [s for s in strings(source) if any(ord(c) > 0xFFFF for c in s)]


def barred_words(source: str) -> list[str]:
    return [s for s in strings(source) if BARRED.search(s)]


def colour_literals(source: str) -> list[str]:
    return [s for s in strings(source) if COLOUR.search(s)]


def _is_tk(node: ast.AST) -> bool:
    return (
        isinstance(node, ast.Attribute)
        and node.attr == "Tk"
        and isinstance(node.value, ast.Name)
        and node.value.id == "tkinter"
    )


def root_class_bindings(source: str) -> list[str]:
    """Module-level names that capture the root class at import, which a patch would miss."""
    found, bound = [], set()
    for stmt in parse(source).body:
        if isinstance(stmt, ast.ImportFrom) and stmt.module == "tkinter":
            for alias in stmt.names:
                if alias.name in ("Tk", "*"):
                    found.append(f"from tkinter import {alias.name}")
                    bound.add(alias.asname or alias.name)
        elif isinstance(stmt, ast.Assign):
            value = stmt.value
            if _is_tk(value) or (isinstance(value, ast.Name) and value.id in bound):
                found.append(ast.unparse(stmt))
                bound.update(t.id for t in stmt.targets if isinstance(t, ast.Name))
    return found


def _imported_names(node: ast.AST) -> list[str]:
    if isinstance(node, ast.Import):
        return [alias.name for alias in node.names]
    if isinstance(node, ast.ImportFrom):
        base = ["eulerase.gui", "eulerase"][node.level - 1] if node.level else ""
        module = ".".join(part for part in (base, node.module) if part)
        return [module, *(f"{module}.{alias.name}" for alias in node.names)]
    return []


def _literal_import(node: ast.AST) -> list[str]:
    if not (isinstance(node, ast.Call) and node.args):
        return []
    func = node.func
    name = func.id if isinstance(func, ast.Name) else getattr(func, "attr", "")
    first = node.args[0]
    if name in DYNAMIC_IMPORTS and isinstance(first, ast.Constant) and isinstance(first.value, str):
        return [first.value]
    return []


def toolkit_imports(source: str) -> list[str]:
    """Every import, caught or not, of a banned module, by statement or by literal call."""
    names = [
        n for node in ast.walk(parse(source)) for n in _imported_names(node) + _literal_import(node)
    ]
    return [n for n in names if any(n == b or n.startswith(f"{b}.") for b in PRESENTER_BANNED)]


def gui_sources() -> dict[Path, str]:
    assert GUI.is_dir(), f"{GUI} is missing"
    sources = {path: path.read_text(encoding="utf-8") for path in sorted(GUI.glob("*.py"))}
    assert GUI / "presenter.py" in sources, "presenter.py is missing"
    return sources


@pytest.mark.parametrize(
    "rule",
    [outside_the_basic_plane, barred_words, colour_literals, root_class_bindings],
    ids=lambda rule: rule.__name__,
)
def test_no_gui_file_breaks_the_rule(rule: Callable[[str], list[str]]):
    broken = {path.name: found for path, source in gui_sources().items() if (found := rule(source))}
    assert not broken


def test_the_presenter_names_no_toolkit_anywhere():
    assert toolkit_imports(gui_sources()[GUI / "presenter.py"]) == []


BAD_SOURCES = [
    (outside_the_basic_plane, 'x = "\U0001f34b"'),
    (barred_words, 'x = "Read Depth"'),
    (barred_words, 'x = f"cov={1}"'),
    (colour_literals, 'x = "#F9D923"'),
    (root_class_bindings, "from tkinter import Tk"),
    (root_class_bindings, "from tkinter import Tk as Root"),
    (root_class_bindings, "from tkinter import *"),
    (root_class_bindings, "import tkinter\nTk = tkinter.Tk"),
    (root_class_bindings, "import tkinter\nTk = tkinter.Tk\nRoot = Tk"),
    (toolkit_imports, "try:\n    import tkinter\nexcept ImportError:\n    pass"),
    (toolkit_imports, "import eulerase.render"),
    (toolkit_imports, "from eulerase import render"),
    (toolkit_imports, "from eulerase.gui import app"),
    (toolkit_imports, "from . import app"),
    (toolkit_imports, "import matplotlib.pyplot"),
    (toolkit_imports, 'import importlib\nimportlib.import_module("tkinter")'),
    (toolkit_imports, '__import__("tkinter")'),
]
GOOD_SOURCES = [
    (outside_the_basic_plane, 'x = "caf\u00e9"'),
    (barred_words, 'x = "times seen"'),
    (colour_literals, 'x = "#fff"'),
    (root_class_bindings, "import tkinter\n\ndef start():\n    return tkinter.Tk()"),
    (
        toolkit_imports,
        "from eulerase import analysis\nimport importlib\nimportlib.import_module(name)",
    ),
]


@pytest.mark.parametrize(
    "rule, source", BAD_SOURCES, ids=lambda x: x if isinstance(x, str) else x.__name__
)
def test_oracle_rejects_a_source_that_breaks_the_rule(rule, source):
    assert rule(source)


@pytest.mark.parametrize(
    "rule, source", GOOD_SOURCES, ids=lambda x: x if isinstance(x, str) else x.__name__
)
def test_oracle_accepts_a_source_that_keeps_the_rule(rule, source):
    assert rule(source) == []


def test_every_rule_fails_when_the_folder_has_no_presenter(monkeypatch, tmp_path):
    monkeypatch.setattr(f"{__name__}.GUI", tmp_path)
    with pytest.raises(AssertionError, match="presenter.py is missing"):
        gui_sources()

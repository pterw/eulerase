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


def _calls_outside(source: str, name: str, home: str | None) -> list[str]:
    """Calls to `name` made anywhere but inside the function `home` (anywhere, if None)."""
    found: list[str] = []

    def visit(node: ast.AST, scope: str) -> None:
        if isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef):
            scope = node.name
        if isinstance(node, ast.Call):
            func = node.func
            if getattr(func, "id", getattr(func, "attr", "")) == name and scope != home:
                found.append(f"{scope}: {ast.unparse(node)}")
        for child in ast.iter_child_nodes(node):
            visit(child, scope)

    visit(parse(source), "<module>")
    return found


def draws_the_figure_itself(source: str) -> list[str]:
    return _calls_outside(source, "draw_figure", None)


def lays_out_outside_relayout(source: str) -> list[str]:
    return _calls_outside(source, "compute_layout", "_relayout")


def builds_a_scene_outside_build_scene(source: str) -> list[str]:
    return _calls_outside(source, "StepScene", "_build_scene")


def _outside_type_checking(node: ast.AST):
    for child in ast.iter_child_nodes(node):
        test = child.test if isinstance(child, ast.If) else None
        if getattr(test, "id", getattr(test, "attr", "")) == "TYPE_CHECKING":
            yield from (n for stmt in child.orelse for n in ast.walk(stmt))
        else:
            yield child
            yield from _outside_type_checking(child)


def runtime_analysis_import(source: str) -> list[str]:
    """Imports of eulerase.analysis that run, which the window gets only through the presenter."""
    names = [n for node in _outside_type_checking(parse(source)) for n in _imported_names(node)]
    return [n for n in names if n == "eulerase.analysis" or n.startswith("eulerase.analysis.")]


ONLY_IN = {runtime_analysis_import: "app.py"}  # the presenter owns the analysis


def gui_sources() -> dict[Path, str]:
    assert GUI.is_dir(), f"{GUI} is missing"
    sources = {path: path.read_text(encoding="utf-8") for path in sorted(GUI.glob("*.py"))}
    assert GUI / "presenter.py" in sources, "presenter.py is missing"
    return sources


@pytest.mark.parametrize(
    "rule",
    [
        outside_the_basic_plane,
        barred_words,
        colour_literals,
        root_class_bindings,
        draws_the_figure_itself,
        lays_out_outside_relayout,
        builds_a_scene_outside_build_scene,
        runtime_analysis_import,
    ],
    ids=lambda rule: rule.__name__,
)
def test_no_gui_file_breaks_the_rule(rule: Callable[[str], list[str]]):
    broken = {
        path.name: found
        for path, source in gui_sources().items()
        if ONLY_IN.get(rule, path.name) == path.name and (found := rule(source))
    }
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
    (draws_the_figure_itself, "def show_step(self):\n    render.draw_figure(self.figure)"),
    (draws_the_figure_itself, "draw_figure(figure)"),
    (lays_out_outside_relayout, "def _build_scene(self):\n    render.compute_layout(a, 1.0)"),
    (builds_a_scene_outside_build_scene, "def show_step(self):\n    render.StepScene(f, a, l)"),
    (
        builds_a_scene_outside_build_scene,
        "def _build_scene(self):\n    def go():\n        StepScene()",
    ),
    (runtime_analysis_import, "from eulerase.analysis import TABLES"),
    (runtime_analysis_import, "from eulerase import analysis"),
    (runtime_analysis_import, "import eulerase.analysis as a"),
    (
        runtime_analysis_import,
        "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    pass\nelse:\n"
        "    from eulerase.analysis import TABLES",
    ),
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
    (
        draws_the_figure_itself,
        "def show_step(self):\n    self.scene.show(1)\n    self.canvas.draw_idle()",
    ),
    (lays_out_outside_relayout, "def _relayout(self):\n    render.compute_layout(a, 1.0)"),
    (builds_a_scene_outside_build_scene, "def _build_scene(self):\n    render.StepScene(f, a, l)"),
    (
        runtime_analysis_import,
        "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from eulerase.analysis import "
        "Analysis\nfrom eulerase import render",
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

"""The command table, the shortcuts and the dialogs of the window, checked without a window."""

import pytest

pytest.importorskip("tkinter", reason="the window module imports tkinter; nothing here opens one")
import tkinter  # noqa: E402

assert getattr(tkinter.Tk, "eulerase_window_guard", False) and getattr(
    tkinter.Toplevel, "eulerase_window_guard", False
), "guard not installed"

from eulerase.gui import app  # noqa: E402
from eulerase.gui.presenter import EXAMPLES, EXPORT_KINDS  # noqa: E402

EDITING_LETTERS = set("acvxzy")


def shortcut_keys(commands) -> list[frozenset[str]]:
    """The key sets of the commands that have a shortcut, case folded."""
    return [
        frozenset(part.lower() for part in c.shortcut.split("+"))
        for c in commands
        if c.shortcut is not None
    ]


def clashes(commands) -> list[frozenset[str]]:
    keys = shortcut_keys(commands)
    return [k for k in keys if keys.count(k) > 1]


def editing_keys(commands) -> list[frozenset[str]]:
    """Ctrl plus a text-editing letter without Shift, or a bare letter or space."""
    found = []
    for keys in shortcut_keys(commands):
        bare = len(keys) == 1 and (next(iter(keys)) == "space" or len(next(iter(keys))) == 1)
        ctrl_letter = (
            len(keys) == 2 and "ctrl" in keys and bool((keys - {"ctrl"}) & EDITING_LETTERS)
        )
        if bare or ctrl_letter:
            found.append(keys)
    return found


def command(shortcut):
    return app.Command("File", "x", shortcut, "export", ("figure",))


def test_shortcuts_are_unique_and_leave_the_text_editing_keys_alone():
    assert sum(c.shortcut is not None for c in app.COMMANDS) > 5
    assert clashes(app.COMMANDS) == []
    assert editing_keys(app.COMMANDS) == []


def test_oracle_rejects_a_clashing_shortcut():
    twins = [command("Ctrl+Shift+F"), command("ctrl+shift+f")]
    assert clashes(twins) != []
    assert editing_keys([command("Ctrl+C")]) != []
    assert editing_keys([command("Ctrl+Shift+C"), command("Escape")]) == []
    assert editing_keys([command("A"), command("Space")]) != []


def test_every_export_and_the_main_actions_have_a_command():
    by_method = {(c.method, c.args): c for c in app.COMMANDS}
    for kind in EXPORT_KINDS:
        assert by_method[("export", (kind,))].shortcut is not None
    for method in (
        "open_file",
        "build",
        "stop",
        "toggle_play",
        "step_back",
        "step_forward",
        "step_start",
        "step_end",
        "quit",
    ):
        assert (method, ()) in by_method
    assert [c.args for c in app.COMMANDS if c.method == "load_example"] == [
        (e.name,) for e in EXAMPLES
    ]


def test_every_command_has_a_label_and_a_menu():
    assert all(c.label and c.menu for c in app.COMMANDS)
    assert {c.menu for c in app.COMMANDS} == {"File", "Examples", "Walk"}


class FakeMenu:
    def __init__(self):
        self.states = {}

    def entryconfigure(self, index, **options):
        self.states[index] = options["state"]


def test_set_exports_enabled_switches_each_entry_by_kind():
    menu = FakeMenu()
    entries = {kind: (menu, i) for i, kind in enumerate(EXPORT_KINDS)}
    wanted = {kind: i % 2 == 0 for i, kind in enumerate(EXPORT_KINDS)}
    assert len(set(wanted.values())) == 2
    window = type("W", (), {"_export_entries": entries})()
    app.Window.set_exports_enabled(window, wanted)
    assert menu.states == {
        i: "normal" if wanted[kind] else "disabled" for i, kind in enumerate(EXPORT_KINDS)
    }


@pytest.mark.parametrize("kind", ["csv", "figure"])
def test_only_the_csv_dialog_skips_the_overwrite_prompt(kind, monkeypatch):
    seen = {}

    def fake_dialog(**options):
        seen.update(options)
        return "chosen"

    monkeypatch.setattr(app.tkinter.filedialog, "asksaveasfilename", fake_dialog)
    window = type("W", (), {"root": object()})()
    assert app.Window.ask_save_path(window, kind, "suggested") == "chosen"
    assert seen["confirmoverwrite"] is (kind != "csv")
    monkeypatch.setattr(app.tkinter.filedialog, "asksaveasfilename", lambda **o: "")
    assert app.Window.ask_save_path(window, kind, "suggested") is None


TK_SEQUENCES = {
    "Ctrl+Shift+O": "<Control-Shift-O>",
    "Ctrl+Shift+F": "<Control-Shift-F>",
    "Ctrl+Shift+C": "<Control-Shift-C>",
    "Ctrl+Shift+A": "<Control-Shift-A>",
    "Ctrl+Shift+B": "<Control-Shift-B>",
    "Ctrl+Shift+G": "<Control-Shift-G>",
    "Ctrl+Q": "<Control-q>",  # without Shift a letter key arrives lower case
    "Ctrl+Enter": "<Control-Return>",
    "Esc": "<Escape>",
    "F5": "<F5>",
    "Alt+Left": "<Alt-Left>",
    "Alt+Right": "<Alt-Right>",
    "Alt+Home": "<Alt-Home>",
    "Alt+End": "<Alt-End>",
}


@pytest.mark.parametrize("shortcut", sorted(TK_SEQUENCES))
def test_each_shortcut_becomes_the_sequence_tk_delivers(shortcut):
    assert shortcut in {c.shortcut for c in app.COMMANDS}
    assert app.tk_sequence(shortcut) == TK_SEQUENCES[shortcut]


def test_every_shortcut_in_the_table_has_a_tested_sequence():
    assert {c.shortcut for c in app.COMMANDS if c.shortcut} == set(TK_SEQUENCES)


@pytest.mark.parametrize("answer", [False, True])
@pytest.mark.parametrize("existing", [[], ["/out/b.csv", "/out/a.csv"]])
def test_confirm_export_always_asks_and_names_the_files_that_exist(existing, answer, monkeypatch):
    asked = []

    def recorder(title, message, **options):
        asked.append(message)
        return answer

    monkeypatch.setattr(app.tkinter.messagebox, "askyesno", recorder)
    paths = ["/out/a.csv", "/out/b.csv"]
    assert existing != paths  # the replaced section is a different list from the written one
    window = type("W", (), {"root": object()})()

    assert app.Window.confirm_export(window, paths, existing) is answer

    (message,) = asked
    written, _, replaced = message.partition("\n\n")
    assert written.splitlines()[1:] == paths
    assert replaced.splitlines()[1:] == existing  # the heading's wording is not pinned

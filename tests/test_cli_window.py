"""The command starts the window for no input or --gui, and only then (never a real window)."""

import subprocess
import sys

import pytest
from window_probe import blocked_code

from eulerase import cli
from eulerase.analysis import DEFAULT_K
from eulerase.errors import EuleraseError

pytest.importorskip("tkinter", reason="the window module imports tkinter; skip where it is absent")

from genomes import child_env  # noqa: E402

from eulerase.gui.presenter import FormValues  # noqa: E402

_STARTED = 7  # neither 0 nor 1, so main is seen to return what start returned


@pytest.fixture
def starts(monkeypatch):
    """Replace the window's start with a fake; the list holds each argument it was given."""
    calls = []

    def fake(form=None):
        calls.append(form)
        return _STARTED

    monkeypatch.setattr("eulerase.gui.app.start", fake, raising=False)
    return calls


def _run(capsys, *argv):
    try:
        code = cli.main(list(argv))
    except SystemExit as stop:
        code = stop.code
    out, err = capsys.readouterr()
    return code, out, err


def test_no_arguments_start_the_window_with_nothing_loaded(starts, capsys):
    assert _run(capsys) == (_STARTED, "", "")
    assert starts == [None]


def _form(text, mode, k, **flags):
    values = {
        "both_strands": False,
        "unique": False,
        "min_multiplicity_text": "",
        "clip_tips": False,
    }
    return FormValues(text, mode, str(k), **{**values, **flags})


@pytest.mark.parametrize(
    ("argv", "expected"),
    [
        (
            ["--seq", "ACGTACGT", "TTGCAACG", "-k", "5", "--both-strands", "--min-count", "2"]
            + ["--clip-tips"],
            _form(
                "ACGTACGT\nTTGCAACG",
                "sequences",
                5,
                both_strands=True,
                min_multiplicity_text="2",
                clip_tips=True,
            ),
        ),
        (
            ["--kmers", "ACG", "CGT", "--clip-tips"],
            _form("ACG CGT", "kmers", DEFAULT_K, clip_tips=True),
        ),
        (
            ["--file", "{fasta}", "--unique"],
            _form(">r1\nACGTACGT\n>r2\nTTGCAACG", "sequences", DEFAULT_K, unique=True),
        ),
    ],
    ids=["sequences-with-settings", "kmers", "file"],
)
def test_gui_with_input_and_settings_loads_them_into_the_window(
    argv, expected, starts, capsys, tmp_path
):
    fasta = tmp_path / "reads.fasta"
    fasta.write_text(">r1\nACGTACGT\n>r2\nTTGCAACG\n", encoding="utf-8")
    assert DEFAULT_K != 5  # the typed k and the default k are told apart
    argv = [str(fasta) if part == "{fasta}" else part for part in argv]
    assert _run(capsys, "--gui", *argv) == (_STARTED, "", "")
    assert starts == [expected]


def test_input_without_gui_never_starts_the_window(starts, capsys):
    code, out, err = _run(capsys, "--seq", "ACGTACGT", "-k", "3")
    assert (code, err) == (0, "")
    assert out.startswith("Eulerian path")
    assert starts == []


@pytest.mark.parametrize(
    ("option", "extra"),
    [
        ("--figure", ["out.png"]),
        ("--csv", ["out"]),
        ("--gfa", ["out.gfa"]),
        ("--contigs", ["out.fa"]),
        ("--genbank", ["out.gb"]),
        ("--metrics", []),
    ],
)
def test_gui_with_an_option_the_window_would_ignore_is_an_error(option, extra, starts, capsys):
    code, out, err = _run(capsys, "--gui", "--seq", "ACGTACGT", option, *extra)
    assert (code, out) == (2, "")
    assert err.startswith("eulerase: ")
    assert len(err.splitlines()) == 1
    assert option in err
    assert starts == []


def test_gui_with_a_setting_and_no_input_is_an_error(starts, capsys):
    code, out, err = _run(capsys, "--gui", "-k", "5")
    assert (code, out) == (2, "")
    assert "-k" in err
    assert starts == []


def test_a_file_that_cannot_be_read_fails_before_the_window_opens(starts, capsys, tmp_path):
    code, out, err = _run(capsys, "--gui", "--file", str(tmp_path / "missing.fasta"))
    assert (code, out) == (1, "")
    assert err.startswith("eulerase: ")
    assert len(err.splitlines()) == 1
    assert starts == []


def _no_toolkit(monkeypatch):
    monkeypatch.setitem(sys.modules, "eulerase.gui.app", None)


def _no_display(monkeypatch):
    def start(form=None):
        raise EuleraseError("The window needs a display.")

    monkeypatch.setattr("eulerase.gui.app.start", start, raising=False)


@pytest.mark.parametrize(
    ("setup", "words"),
    [(_no_toolkit, ["tkinter", "still work"]), (_no_display, ["needs a display"])],
    ids=["no-toolkit", "no-display"],
)
def test_a_missing_toolkit_or_display_is_one_line_and_exit_1(setup, words, monkeypatch, capsys):
    setup(monkeypatch)
    code, out, err = _run(capsys)
    assert (code, out) == (1, "")
    assert err.startswith("eulerase: ")
    assert len(err.splitlines()) == 1
    assert "Traceback" not in err
    assert all(word in err for word in words)


def test_running_with_input_does_not_import_the_toolkit():
    # blocked_code makes any import of tkinter raise, so a mistaken import of the window
    # fails the child, and no window can open. The block stores None, so a plain `in` is true.
    body = (
        "from eulerase.cli import main\n"
        "assert main(['--seq', 'ACGTACGT', '-k', '3']) == 0\n"
        "assert sys.modules.get('tkinter') is None\n"
    )
    env = child_env()
    env["MPLBACKEND"] = "Agg"
    done = subprocess.run(
        [sys.executable, "-c", blocked_code(body)],
        env=env,
        capture_output=True,
        timeout=120,
        check=False,
    )
    assert done.returncode == 0, done.stderr.decode()

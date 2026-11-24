"""A file argument that names nothing readable is an error, never a sequence."""

import os

import pytest

from eulerase import cli


def _fails_cleanly(code: int, out: str, err: str) -> bool:
    """Exit 1, one line on stderr naming the problem, nothing on standard output."""
    lines = err.splitlines()
    return code == 1 and out == "" and len(lines) == 1 and lines[0].startswith("eulerase: ")


def _run_file(capsys, path: str) -> tuple[int, str, str]:
    code = cli.main(["--file", path, "-k", "3"])
    out, err = capsys.readouterr()
    return code, out, err


@pytest.mark.parametrize("kind", ["plain name", "missing folder", "sequence text"])
def test_missing_file_is_an_error(kind, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    path = {
        "plain name": "no_such.fa",
        "missing folder": str(tmp_path / "gone" / "reads.fa"),
        # Valid DNA at k=3: were it read as a sequence, the run would succeed.
        "sequence text": "ACGTACGT",
    }[kind]
    assert not os.path.exists(path)
    assert _fails_cleanly(*_run_file(capsys, path))


def test_directory_given_as_file_is_an_error(tmp_path, capsys):
    assert tmp_path.is_dir()
    assert _fails_cleanly(*_run_file(capsys, str(tmp_path)))


def test_oracle_rejects_a_file_read_as_a_sequence(capsys):
    # What a path taken as a sequence would print: a result, exit 0, nothing on stderr.
    assert cli.main(["--seq", "ACGTACGT", "-k", "3"]) == 0
    out, err = capsys.readouterr()
    assert out != "" and err == ""
    assert not _fails_cleanly(0, out, err)
    assert not _fails_cleanly(1, "", "eulerase: a\neulerase: b\n")
    assert not _fails_cleanly(1, "a result", "eulerase: a\n")
    assert _fails_cleanly(1, "", "eulerase: a\n")

"""The eulerase command: the input forms, the results printed, and how failures are reported."""

import subprocess
import sys

import pytest
from genomes import brute_kmers, child_env

from eulerase import cli
from eulerase.analysis import DEFAULT_K

_RUN_MAIN = "import sys; from eulerase.cli import main; sys.exit(main(sys.argv[1:]))"


def _run(capsys, *argv: str) -> tuple[int, str, str]:
    code = cli.main(list(argv))
    out, err = capsys.readouterr()
    return code, out, err


def _child(*argv: str) -> subprocess.CompletedProcess[bytes]:
    """Run main in a fresh interpreter whose standard output starts as cp1252."""
    env = child_env()
    env["PYTHONIOENCODING"] = "cp1252"
    return subprocess.run(
        [sys.executable, "-c", _RUN_MAIN, *argv],
        env=env,
        capture_output=True,
        timeout=120,
        check=False,
    )


def test_prints_results_for_the_ten_copies_case():
    # ATGC ten times at k=3: k-mers ATG and TGC seen ten times each, so 20 edges over the
    # nodes AT, TG, GC; ten trails of two edges, each spelling ATGC.
    expected = (
        "No single Eulerian path · minimum cover = 10 trails · 20 edges · 3 nodes\n"
        + "".join(f"trail {i} (2 edges): ATGC\n" for i in range(1, 11))
    )
    done = _child("--seq", *["ATGC"] * 10)
    assert done.returncode == 0, done.stderr
    assert b"\xc2\xb7" in expected.encode("utf-8")  # the middle dot is two bytes in UTF-8
    # Only the cp1252 half can fail on Linux and macOS; the \r half can only fail on Windows.
    assert done.stdout == expected.encode("utf-8")
    assert done.stderr == b""


def test_main_prints_non_ascii_input_as_utf8():
    # Six Greek letters, long enough for the default k, none of them writable in cp1252, so
    # a stream left as cp1252 would fail. A record name is never printed, so only the
    # sequence can show the output encoding. The parser upper-cases the letters.
    word = "αβγδεζ"
    with pytest.raises(UnicodeEncodeError):
        word.encode("cp1252")
    done = _child("--seq", word)
    assert done.returncode == 0, done.stderr
    assert f"trail 1 (4 edges): {word.upper()}\n".encode() in done.stdout


def test_main_writes_an_error_sentence_as_utf8():
    # A Greek file name that does not exist, none of its letters writable in cp1252, so a
    # standard error left as cp1252 would escape them instead of writing the name.
    name = "αβγ.fa"
    with pytest.raises(UnicodeEncodeError):
        name.encode("cp1252")
    done = _child("--file", name)
    assert done.returncode == 1
    assert done.stderr == f"eulerase: The file '{name}' does not exist; check the path.\n".encode()


def test_each_input_form_gives_the_same_trails(tmp_path, capsys):
    # No two-letter pair repeats, so the one Eulerian path spells the input itself.
    seq = "AACCGGTTAGC"
    pairs = [seq[i : i + 2] for i in range(len(seq) - 1)]
    assert len(set(pairs)) == len(pairs)
    fasta = tmp_path / "reads.fasta"
    fasta.write_text(f">one\n{seq[:6]}\n{seq[6:]}\n", encoding="utf-8")
    expected = f"Eulerian path · 9 edges · 10 nodes\ntrail 1 (9 edges): {seq}\n"
    assert brute_kmers([seq], 3, False)[0] == "AAC"  # nine k-mers at the default k
    for argv in (
        ["--seq", seq],
        ["--kmers", *brute_kmers([seq], 3, False)],
        ["--file", str(fasta)],
    ):
        code, out, err = _run(capsys, *argv)
        assert (code, out, err) == (0, expected, "")


@pytest.mark.parametrize(
    "argv",
    [
        ["--seq", "ACGT", "-k", "1"],
        ["--seq", ""],
        ["--seq", "HELLO", "--both-strands"],
        ["--file", ""],
        ["--file", "no_such.fa", "-k", "3"],
        ["--file", "gone/reads.fa", "-k", "3"],
        # Valid DNA at k=3: were it read as a sequence, the run would succeed.
        ["--file", "ACGTACGT", "-k", "3"],
        ["--file", ".", "-k", "3"],
    ],
    ids=[
        "k-too-small",
        "empty-sequence",
        "both-strands-not-dna",
        "empty-file-path",
        "missing-file",
        "missing-folder",
        "sequence-text",
        "directory",
    ],
)
def test_a_failure_is_one_line_on_stderr_and_exit_1(argv, capsys, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # a folder holding none of the files and folders named above
    code, out, err = _run(capsys, *argv)
    assert code == 1
    assert out == ""
    assert err.startswith("eulerase: ")
    assert err.endswith("\n")
    assert len(err.splitlines()) == 1
    assert "Traceback" not in err


@pytest.mark.parametrize(
    "argv", [["--bogus"], ["--seq", "ACGT", "-k"], ["--seq", "ACGT", "-k", "three"]]
)
def test_usage_errors_are_one_line_and_exit_2(argv, capsys):
    with pytest.raises(SystemExit) as stop:
        cli.main(argv)
    out, err = capsys.readouterr()
    assert stop.value.code == 2
    assert out == ""
    assert err.startswith("eulerase: ")
    assert len(err.splitlines()) == 1


def test_k_defaults_to_the_shared_default(capsys):
    # At k=3 and k=4 this sequence gives different trail counts, so a wrong default shows.
    seq = "ATGCATGCATGC"
    _, default_out, _ = _run(capsys, "--seq", seq)
    _, given_out, _ = _run(capsys, "--seq", seq, "-k", str(DEFAULT_K))
    _, other_out, _ = _run(capsys, "--seq", seq, "-k", str(DEFAULT_K + 1))
    assert other_out != given_out
    assert default_out == given_out


def test_an_error_that_is_not_a_eulerase_error_is_not_hidden(monkeypatch):
    def broken(records, settings):
        raise RuntimeError("a bug, not a user mistake")

    monkeypatch.setattr(cli, "analyze", broken)
    with pytest.raises(RuntimeError):
        cli.main(["--seq", "ACGTACGT"])

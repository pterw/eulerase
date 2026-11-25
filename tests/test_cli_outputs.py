"""The files the command writes: each option, its bytes, the CSV set and the order of the lines."""

import csv
import io
import re
import sys
from pathlib import Path

import matplotlib
import pytest
from Bio import SeqIO
from genomes import random_dna

from eulerase import cli
from eulerase import io as eio
from eulerase.analysis import Settings, analyze
from eulerase.sequence import parse_records

_SEQ = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"  # 39 letters, ORF-bearing DNA
_PNG = b"\x89PNG\r\n\x1a\n"


def _analysis():
    return analyze(parse_records(_SEQ), Settings(k=4))


def _run(capsys, *argv: str) -> list[str]:
    assert cli.main(list(argv)) == 0
    out, err = capsys.readouterr()
    assert err == ""
    return out.splitlines()


def _paths(tmp_path: Path) -> dict[str, str]:
    return {
        "figure": str(tmp_path / "graph.png"),
        "csv": str(tmp_path / "tables"),
        "gfa": str(tmp_path / "out.gfa"),
        "contigs": str(tmp_path / "out.fasta"),
        "genbank": str(tmp_path / "out.gb"),
    }


def _all_options(paths: dict[str, str]) -> list[str]:
    return [
        "--seq", _SEQ, "-k", "4",
        "--figure", paths["figure"],
        "--csv", paths["csv"],
        "--gfa", paths["gfa"],
        "--contigs", paths["contigs"],
        "--genbank", paths["genbank"],
    ]  # fmt: skip


def test_each_output_option_writes_its_file(tmp_path, capsys):
    paths = _paths(tmp_path)
    _run(capsys, *_all_options(paths))
    a = _analysis()

    fasta = list(SeqIO.parse(paths["contigs"], "fasta"))
    assert [(r.id, str(r.seq)) for r in fasta] == [(f"ctg{c.id}", c.sequence) for c in a.contigs]
    assert len(fasta) >= 1

    genbank = list(SeqIO.parse(paths["genbank"], "genbank"))
    assert [str(r.seq) for r in genbank] == [c.sequence for c in a.contigs]

    segments = [
        line.split("\t") for line in Path(paths["gfa"]).read_text("utf-8").splitlines()
        if line.startswith("S\t")
    ]  # fmt: skip
    assert [(s[1], s[2]) for s in segments] == [(f"ctg{c.id}", c.sequence) for c in a.contigs]

    edges = list(csv.reader(io.StringIO(Path(f"{paths['csv']}_edges.csv").read_text("utf-8"))))
    assert edges[0][3] == "k-mer"
    assert sorted(row[3] for row in edges[1:]) == sorted(
        _SEQ[i : i + 4] for i in range(len(_SEQ) - 3)
    )

    figure = Path(paths["figure"]).read_bytes()
    assert figure.startswith(_PNG)
    assert len(figure) > 5000


def test_written_files_are_byte_identical_to_the_writers(tmp_path, capsys):
    paths = _paths(tmp_path)
    _run(capsys, *_all_options(paths))
    a = _analysis()
    expected = {
        paths["contigs"]: eio.fasta_text(a),
        paths["genbank"]: eio.genbank_text(a),
        paths["gfa"]: eio.gfa_text(a),
    }
    for path, text in expected.items():
        data = Path(path).read_bytes()
        assert data == text.encode("utf-8")
        assert b"\r" not in data
        assert not data.startswith(b"\xef\xbb\xbf")
        assert text  # an empty file would match an empty writer


def test_csv_writes_the_five_tables(tmp_path, capsys):
    paths = _paths(tmp_path)
    _run(capsys, *_all_options(paths))
    a = _analysis()
    written = eio.csv_paths(paths["csv"])
    assert sorted(written) == ["assembly", "contigs", "edges", "nodes", "trails"]
    for table, path in written.items():
        assert Path(path).read_bytes() == eio.csv_text(a, table.capitalize()).encode("utf-8")
    assert len({Path(p).read_bytes() for p in written.values()}) == 5  # each table is its own


def test_wrote_lines_come_in_the_frozen_order_and_print_each_path_as_typed(
    tmp_path, capsys, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "sub").mkdir()
    # Relative, so a line that resolved or made the path absolute would differ from what was typed.
    paths = {
        "figure": "graph.png",
        "csv": "./tables",
        "gfa": "out.gfa",
        "contigs": "sub/out.fasta",
        "genbank": "sub/out.gb",
    }
    # Typed in an order unlike the output order, so the order cannot come from the command line.
    argv = [
        "--seq", _SEQ, "-k", "4",
        "--genbank", paths["genbank"],
        "--contigs", paths["contigs"],
        "--gfa", paths["gfa"],
        "--csv", paths["csv"],
        "--figure", paths["figure"],
    ]  # fmt: skip
    lines = _run(capsys, *argv)
    csv_files = list(eio.csv_paths(paths["csv"]).values())
    expected = [paths["figure"], *csv_files, paths["gfa"], paths["contigs"], paths["genbank"]]
    assert len(set(expected)) == len(expected) == 9
    assert lines[-9:] == [f"Wrote {p}" for p in expected]
    assert not any(line.startswith("Wrote ") for line in lines[:-9])


def test_figure_for_a_graph_over_the_drawing_limit_says_why(tmp_path, capsys):
    over = tmp_path / "big.svg"
    sequence = random_dna(700, 7)
    nodes = len(analyze(parse_records(sequence), Settings(k=9)).nodes)
    assert nodes > 400
    with matplotlib.rc_context({"svg.fonttype": "none"}):  # keep the text readable in the file
        lines = _run(capsys, "--seq", sequence, "-k", "9", "--figure", str(over))
    assert lines[-1] == f"Wrote {over}"
    svg = over.read_text("utf-8")
    # Only the drawn text counts: a full drawing has "400" in coordinates and the node count in
    # its headline, but no one text element that holds both.
    texts = re.findall(r"<text[^>]*>([^<]*)</text>", svg)
    assert any(f"{nodes}" in t and "400" in t for t in texts)


def test_failed_output_is_an_error_after_the_earlier_files(tmp_path, capsys):
    """A GenBank file cannot hold protein letters: exit 1, one line, the FASTA already written."""
    fasta = tmp_path / "out.fasta"
    code = cli.main(
        ["--seq", "MKVLAAGIV", "-k", "3", "--contigs", str(fasta), "--genbank", str(tmp_path / "g")]
    )
    out, err = capsys.readouterr()
    assert code == 1
    assert err.startswith("eulerase: ") and err.count("\n") == 1
    assert out.splitlines()[-1] == f"Wrote {fasta}"
    assert fasta.exists() and not (tmp_path / "g").exists()


def test_oracle_rejects_a_wrong_file_order(tmp_path, capsys):
    paths = _paths(tmp_path)
    lines = _run(capsys, *_all_options(paths))
    csv_files = list(eio.csv_paths(paths["csv"]).values())
    wrong = [paths["gfa"], paths["figure"], *csv_files, paths["contigs"], paths["genbank"]]
    assert lines[-9:] != [f"Wrote {p}" for p in wrong]


_BAD_OPTIONS = ("--figure", "--csv", "--gfa", "--contigs", "--genbank")
_BAD_NAMES = {
    "--figure": "x.png",
    "--csv": "x",
    "--gfa": "x.gfa",
    "--contigs": "x.fasta",
    "--genbank": "x.gb",
}


def _assert_refused_before_any_work(tmp_path, capsys, monkeypatch, option, run, bad=None):
    """`run(argv)` must refuse a bad `option` path with one stderr line, no output, no file.

    `bad` replaces the missing-folder path with its own arguments and the exact error line.
    """
    calls = []

    def spy(*args, **kwargs):
        calls.append(args)
        raise AssertionError("the analysis ran before the output paths were checked")

    monkeypatch.setattr(cli, "analyze", spy)
    missing = tmp_path / "missing"
    arguments, line = bad or ([option, str(missing / _BAD_NAMES[option])], None)
    # The good option is a different one that the command writes first, where the order allows.
    good_option = "--csv" if option == "--figure" else "--figure"
    good = str(tmp_path / ("good" if good_option == "--csv" else "good.png"))
    good_file = Path(eio.csv_paths(good)["edges"] if good_option == "--csv" else good)
    assert not missing.exists() and not good_file.exists() and good not in arguments

    code = run(["--seq", _SEQ, "-k", "4", good_option, good, *arguments])
    out, err = capsys.readouterr()

    assert code == 1
    assert err.startswith("eulerase: ") and err.count("\n") == 1 and option in err
    assert line is None or err == line
    assert out == ""
    assert not good_file.exists()
    assert calls == []


@pytest.mark.parametrize(
    "option", [*_BAD_OPTIONS, "same file", *(f"{o} empty" for o in _BAD_OPTIONS)]
)
def test_a_bad_output_path_fails_before_any_output_is_written(
    tmp_path, capsys, monkeypatch, option
):
    bad = None
    empty = option.endswith(" empty")
    if empty:
        # An empty path is falsy: a test by truth would skip it and write into the current folder.
        option = option.removesuffix(" empty")
        line = f"eulerase: {option} needs a file name, not an empty path; give one.\n"
        bad = ([option, ""], line)
        monkeypatch.chdir(tmp_path)
    elif option == "same file":
        # The second spelling reaches the first csv-set file through "..": same file, other text.
        clash = f"{tmp_path}/elsewhere/../tables_edges.csv"
        line = (
            f"eulerase: --csv and --genbank name the same file {clash}; "
            "give each output its own path.\n"
        )
        option, bad = "--csv", (["--csv", str(tmp_path / "tables"), "--genbank", clash], line)
    _assert_refused_before_any_work(tmp_path, capsys, monkeypatch, option, cli.main, bad)
    if empty:
        assert list(tmp_path.iterdir()) == []  # nothing went into the current folder


def test_a_folder_that_cannot_be_written_is_refused(tmp_path, capsys, monkeypatch):
    monkeypatch.setattr(cli.os, "access", lambda *a, **k: False)
    code = cli.main(["--seq", _SEQ, "-k", "4", "--gfa", str(tmp_path / "a.gfa")])
    out, err = capsys.readouterr()
    assert code == 1 and out == "" and err.count("\n") == 1 and "--gfa" in err
    assert not (tmp_path / "a.gfa").exists()


def test_an_existing_file_that_cannot_be_written_is_refused(tmp_path, capsys, monkeypatch):
    existing = tmp_path / "out.gfa"
    existing.write_bytes(b"old bytes")
    # Only the file is refused: the folder, which is writable, cannot be what the command meets.
    monkeypatch.setattr(cli.os, "access", lambda path, mode: Path(path) != existing)
    assert cli.os.access(tmp_path, cli.os.W_OK)
    assert not cli.os.access(existing, cli.os.W_OK)

    code = cli.main(["--seq", _SEQ, "-k", "4", "--gfa", str(existing)])
    out, err = capsys.readouterr()
    assert code == 1 and out == "" and err.count("\n") == 1 and "--gfa" in err
    assert existing.read_bytes() == b"old bytes"


@pytest.mark.parametrize("name", ["graph.bmp", "graph"])
def test_unsupported_figure_format_is_refused(tmp_path, capsys, name):
    target = tmp_path / name
    code = cli.main(
        ["--seq", _SEQ, "-k", "4", "--gfa", str(tmp_path / "a.gfa"), "--figure", str(target)]
    )
    out, err = capsys.readouterr()
    assert code == 1
    assert err.startswith("eulerase: ") and err.count("\n") == 1 and "--figure" in err
    assert out == ""
    assert not target.exists() and not (tmp_path / "a.gfa").exists()


def test_output_path_that_is_a_directory_is_refused(tmp_path, capsys):
    folder = tmp_path / "out.gfa"
    folder.mkdir()
    code = cli.main(
        ["--seq", _SEQ, "-k", "4", "--figure", str(tmp_path / "g.png"), "--gfa", str(folder)]
    )
    out, err = capsys.readouterr()
    assert code == 1
    assert err.startswith("eulerase: ") and err.count("\n") == 1 and "--gfa" in err
    assert out == "" and not (tmp_path / "g.png").exists()


def test_a_write_that_fails_after_the_preflight_is_one_line(tmp_path, capsys, monkeypatch):
    def refuse(path, text):
        raise PermissionError("access denied")

    monkeypatch.setattr(eio, "write_text_file", refuse)
    target = str(tmp_path / "out.gfa")
    code = cli.main(["--seq", _SEQ, "-k", "4", "--gfa", target])
    out, err = capsys.readouterr()
    assert code == 1
    assert err == f"eulerase: cannot write {target}: access denied\n"
    assert "Traceback" not in err


def test_oracle_rejects_a_late_failure(tmp_path, capsys, monkeypatch):
    def writes_then_fails(argv):
        Path(argv[argv.index("--figure") + 1]).write_bytes(b"x")
        print("eulerase: late", file=sys.stderr)
        return 1

    def analyses_then_fails(argv):
        cli.analyze([], None)  # the spy that the check installs
        return 1

    for stub in (writes_then_fails, analyses_then_fails):
        with pytest.raises(AssertionError):
            _assert_refused_before_any_work(tmp_path, capsys, monkeypatch, "--gfa", stub)

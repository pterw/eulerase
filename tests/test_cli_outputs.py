"""The files the command writes: each option, its bytes, the CSV set and the order of the lines."""

import csv
import io
from pathlib import Path

import matplotlib
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


def test_wrote_lines_come_in_the_frozen_order(tmp_path, capsys):
    paths = _paths(tmp_path)
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
    text = over.read_text("utf-8")
    assert f"{nodes}" in text and "400" in text


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

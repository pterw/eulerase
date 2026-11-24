"""Options that would be ignored are errors; options that apply still run."""

import pytest

from eulerase import cli

_NEEDS_INPUT = [
    (["-k", "5"], "-k"),
    (["--both-strands"], "--both-strands"),
    (["--unique"], "--unique"),
    (["--min-count", "2"], "--min-count"),
    (["--min-count", "0"], "--min-count"),  # zero is a value the user gave
    (["--clip-tips"], "--clip-tips"),
    (["--metrics"], "--metrics"),
    (["--figure", "out.png"], "--figure"),
    (["--csv", "out"], "--csv"),
    (["--gfa", "out.gfa"], "--gfa"),
    (["--contigs", "out.fa"], "--contigs"),
    (["--genbank", "out.gb"], "--genbank"),
]
_REFUSED = [
    *_NEEDS_INPUT,
    (["--seq", "A", "--kmers", "AC"], "--kmers"),
    (["--seq", "A", "--file", "f"], "--file"),
    (["--kmers", "AC", "AG", "-k", "3"], "-k"),
]


def _assert_refused(main, argv, named, capsys, tmp_path):
    """The command exits 2 with one stderr line that names the option, and does nothing."""
    try:
        code = main(argv)
    except SystemExit as stop:
        code = stop.code
    out, err = capsys.readouterr()
    assert code == 2
    assert out == ""
    assert err.startswith("eulerase: ")
    assert len(err.splitlines()) == 1
    assert named in err
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(("argv", "named"), _REFUSED, ids=[" ".join(a) for a, _ in _REFUSED])
def test_option_that_would_be_ignored_is_an_error(argv, named, capsys, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _assert_refused(cli.main, argv, named, capsys, tmp_path)


def test_every_parser_option_is_in_the_table_of_refused_options():
    options = {s for a in cli._parser()._actions for s in a.option_strings if s.startswith("-")}
    inputs = {"--seq", "--kmers", "--file", "--gui", "-h", "--help"}
    assert options - inputs == {named for _, named in _NEEDS_INPUT}


def test_valid_combinations_still_run(capsys, tmp_path):
    seq = "AACCGGTTAGC"
    fasta = tmp_path / "reads.fasta"
    fasta.write_text(f">one\n{seq}\n", encoding="utf-8")
    # -k 4 gives len(seq) - 3 = 8 edges; the default k of 3 would give 9, so -k is seen to act.
    assert cli.DEFAULT_K == 3
    runs = [
        (
            f"--seq {seq} -k 4 --both-strands --unique --min-count 1 --clip-tips --metrics".split(),
            "Eulerian",
        ),
        (["--kmers", "ACG", "CGT", "--both-strands"], "Eulerian"),
        (["--file", str(fasta), "-k", "4"], "Eulerian path · 8 edges"),
    ]
    for argv, headline in runs:
        assert cli.main(argv) == 0
        out, err = capsys.readouterr()
        assert err == ""
        assert headline in out.splitlines()[0]


def test_oracle_rejects_a_silent_ignore(capsys, tmp_path):
    def ignoring_main(argv):
        """A launcher that drops -k and its value when --kmers is given, as the source did."""
        at = argv.index("-k")
        return cli.main(argv[:at] + argv[at + 2 :])

    with pytest.raises(AssertionError):
        _assert_refused(ignoring_main, ["--kmers", "AC", "AG", "-k", "3"], "-k", capsys, tmp_path)

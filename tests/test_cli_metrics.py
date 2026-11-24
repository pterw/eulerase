"""The notes and the metrics the command prints, and what its help says."""

import re
from collections import Counter

import pytest
from genomes import brute_kmers

from eulerase import cli
from eulerase.analysis import Settings, analyze, rows
from eulerase.sequence import parse_records

_CYCLE = "ATGCATGCATGC"
_LINE = "AACCGGTTAGC"
_CLIP_SENTENCE = "It removes dead ends, not bubbles: for reads, also give --min-count."
_BARRED = re.compile(r"coverage|cov=|depth|bubble", re.IGNORECASE)


def _run(capsys, *argv: str) -> list[str]:
    assert cli.main(list(argv)) == 0
    out, err = capsys.readouterr()
    assert err == ""
    return out.splitlines()


def _measures(lines: list[str]) -> dict[str, str]:
    """The `measure: value` lines of a --metrics run, read back by name."""
    skipped = ("Note: ", "trail ", "contig ")
    return dict(
        line.split(": ", 1) for line in lines if ": " in line and not line.startswith(skipped)
    )


def _flat(text: str) -> str:
    """White space collapsed (argparse wraps --help) and the one allowed denial taken out."""
    return " ".join(text.split()).replace("not bubbles", "")


def _helps() -> list[str]:
    return [a.help for a in cli._parser()._actions if a.help]


def test_metrics_block_matches_the_assembly_table(capsys):
    table = rows(analyze(parse_records(_CYCLE), Settings()), "Assembly")
    expected = [f"{measure}: {value}" for measure, value in table]
    expected.append("contig 1 length 4 multiplicity 2.5000 cycle yes")
    lines = _run(capsys, "--seq", _CYCLE, "--metrics")
    assert lines[-len(expected) :] == expected
    line_lines = _run(capsys, "--seq", _LINE, "--metrics")
    assert line_lines[-1] == "contig 1 length 11 multiplicity 1.0000 cycle no"


@pytest.mark.parametrize(
    ("flags", "measure", "expected"),
    [
        (["-k", "4"], "total k-mers", len(brute_kmers([_CYCLE], 4, False))),
        (["--both-strands"], "total k-mers", len(brute_kmers([_CYCLE], 3, True))),
        (["--unique"], "total k-mers", len(set(brute_kmers([_CYCLE], 3, False)))),
        (
            ["--min-count", "3"],
            "removed k-mers",
            sum(n < 3 for n in Counter(brute_kmers([_CYCLE], 3, False)).values()),
        ),
    ],
    ids=["k", "both-strands", "unique", "min-count"],
)
def test_settings_reach_the_analysis(capsys, flags, measure, expected):
    # Each independent count differs from the default run's, so a flag that is ignored shows.
    default = _measures(_run(capsys, "--seq", _CYCLE, "--metrics"))[measure]
    given = _measures(_run(capsys, "--seq", _CYCLE, *flags, "--metrics"))[measure]
    assert str(expected) != default
    assert given == str(expected)


def test_clip_tips_removes_exactly_the_k_mers_the_contigs_lose(capsys):
    # Three copies of one read and one read that leaves it at a branch: a dead end the clipping
    # removes. The contigs are all linear, so each holds length - (k - 1) distinct k-mers.
    reads = ["CCATGGTCAACG"] * 3 + ["CCATGGTCTT"]
    distinct = len(set(brute_kmers(reads, 3, False)))

    def run(*flags: str) -> list[str]:
        return _run(capsys, "--seq", *reads, *flags, "--metrics")

    plain, clipped = run(), run("--clip-tips")
    assert _measures(plain)["clipped k-mers"] == "0"
    removed = int(_measures(clipped)["clipped k-mers"])
    assert removed > 0
    contig_lines = [line.split() for line in clipped if line.startswith("contig ")]
    assert all(parts[-1] == "no" for parts in contig_lines)
    held = sum(int(parts[3]) - 2 for parts in contig_lines)
    assert held == distinct - removed


def test_notes_are_printed_one_per_line_under_their_condition(capsys):
    # Both notes' conditions are independent facts about the input: a record shorter than k, and
    # a cycle.
    with_cycle = _run(capsys, "--seq", _CYCLE, "AC")
    without = _run(capsys, "--seq", _LINE)
    assert sum(line.startswith("Note: ") for line in with_cycle) >= 2
    assert not any(line.startswith("Note: ") for line in without)


def test_no_output_text_or_help_names_a_count_as_coverage_or_depth(capsys):
    outputs = [
        "\n".join(_run(capsys, "--seq", _CYCLE, *flags, "--metrics"))
        for flags in (
            ["--both-strands", "--min-count", "2", "--clip-tips"],
            ["--unique", "--min-count", "1", "--clip-tips"],
        )
    ]
    with pytest.raises(SystemExit):
        cli.main(["--help"])
    outputs.append(capsys.readouterr().out)
    outputs.extend(_helps())
    assert all(_flat(text) for text in outputs)
    assert not [text for text in outputs if _BARRED.search(_flat(text))]


def test_oracle_rejects_a_help_that_claims_bubble_removal():
    assert _BARRED.search(_flat("It removes bubbles."))
    assert _BARRED.search(_flat("It removes dead ends, not bubbles, and bubbles."))
    assert not _BARRED.search(_flat(_CLIP_SENTENCE))


def test_clip_tips_help_says_what_it_does_not_remove():
    clip = next(a.help for a in cli._parser()._actions if "--clip-tips" in a.option_strings)
    assert " ".join(clip.split()).count(_CLIP_SENTENCE) == 1

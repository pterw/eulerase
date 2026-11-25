"""Sequence helpers: parsing, alphabet check, reverse complement, GC content."""

import random

import pytest
from genomes import revcomp

from eulerase.assembly import assembly_metrics, build_assembly_graph, extract_contigs
from eulerase.kmers import count_kmers
from eulerase.sequence import (
    Record,
    gc_content,
    is_dna,
    parse_kmer_list,
    parse_records,
    reverse_complement,
)


def _random_text(seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGTNacgtn") for _ in range(rng.randint(1, 60)))


@pytest.mark.parametrize("seed", range(20))
def test_reverse_complement_matches_independent_translation(seed):
    seq = _random_text(seed)
    assert reverse_complement(seq).upper() == revcomp(seq.upper())


def test_oracle_rejects_a_wrong_reverse_complement():
    # The independent oracle must tell the reverse from the reverse complement.
    assert revcomp("AACG") == "CGTT"
    assert revcomp("AACG") != "GCAA"
    assert revcomp("AACGN") == "NCGTT"


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        (">a\nACGT\nTTGG\n>b\nCC\n", [Record("a", "ACGTTTGG"), Record("b", "CC")]),
        (">a\r\nacgt\r\n\r\ntt\r\n", [Record("a", "ACGTTT")]),
        ("\n\n>a\n\nAC\n\n\nGT\n\n", [Record("a", "ACGT")]),
        ("acgt\n\nTTGG\r\n", [Record("seq1", "ACGT"), Record("seq2", "TTGG")]),
        ("AC GT\n", [Record("seq1", "ACGT")]),
        (">only a header\n", [Record("only a header", "")]),  # was [] before this was inverted
        (">a\n>b\nAC\n", [Record("a", ""), Record("b", "AC")]),
        ("", []),
    ],
)
def test_parse_records_normalises_every_layout(text, expected):
    assert parse_records(text) == expected


@pytest.mark.parametrize(
    "text",
    [
        ">ACG\nAAC, GGT;TTA\nCCA",
        ">ACG\r\naac,ggt;tta\r\nCCA\r\n",
        "  AAC GGT\t TTA ,, ;CCA ",
    ],
)
def test_parse_kmer_list_drops_header_lines_and_splits(text):
    # The header word ACG has the length of a k-mer and must never become one.
    names = [r.name for r in parse_kmer_list(text)]
    assert [r.sequence for r in parse_kmer_list(text)] == ["AAC", "GGT", "TTA", "CCA"]
    assert names == ["kmer1", "kmer2", "kmer3", "kmer4"]


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [("ACGTN", True), ("acgt", True), ("LEMON", False), ("", False), ("ACGU", False)],
)
def test_is_dna_accepts_acgtn_only(sequence, expected):
    assert is_dna(sequence) is expected


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [
        ("ACGT", 0.5),
        ("GCAT", 0.5),
        ("gcat", 0.5),
        ("GGCC", 1.0),
        ("GCNN", 1.0),
        ("ACG", 2 / 3),
        ("GCNNA", 2 / 3),
        ("NNNN", None),
        ("", None),
        ("LEMON", None),
        ("ACGU", None),
    ],
)
def test_gc_content_cases(sequence, expected):
    assert gc_content(sequence) == expected


def _hand_gc(sequence: str) -> float | None:
    # Independent count of (G + C) / (A + C + G + T); None for non-DNA or no known base.
    upper = sequence.upper()
    if not upper or set(upper) - set("ACGTN"):
        return None
    known = sum(upper.count(base) for base in "ACGT")
    return (upper.count("G") + upper.count("C")) / known if known else None


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [("GCNNA", 2 / 3), ("ACGTGC", 4 / 6), ("LEMONLEMONADE", None)],
)
def test_every_caller_gets_the_same_gc(sequence, expected):
    counts = count_kmers([sequence[i : i + 3] for i in range(len(sequence) - 2)])
    contigs = extract_contigs(build_assembly_graph(counts), 3, both_strands=False)
    assert contigs
    for contig in contigs:
        assert contig.gc == assembly_metrics([contig]).gc == gc_content(contig.sequence)
        assert contig.gc == _hand_gc(contig.sequence)
    if len(contigs) == 1:
        assert contigs[0].sequence == sequence and contigs[0].gc == expected
    else:
        assert expected is None

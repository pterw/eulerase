"""No output calls an occurrence count "coverage" or "depth" (a check for barred words)."""

import pytest
from genomes import random_dna

from eulerase.analysis import TABLES, Settings, analyze, columns
from eulerase.io import csv_text, fasta_text, genbank_text, gfa_text
from eulerase.sequence import Record

BARRED = ("coverage", "cov=", "depth")


def barred_in(text: str) -> list[str]:
    return [word for word in BARRED if word in text.lower()]


_MAIN = "AGGAGTTAAATCGATGTCTCCTTCTGGCTT"
_BRANCH = "AGGAGTTAAATCGATGCCA"

# (reads, settings, the kinds of note the analysis must give: one note each)
_CASES = (
    # GCA and CAT are seen twice, so a minimum multiplicity of 3 removes them (plural form)
    (["ATGCATGCATGC"], Settings(k=3, min_multiplicity=3), ("minimum",)),
    (["ACGT", "ACGT", "ACGA"], Settings(k=3, min_multiplicity=2), ("minimum, one",)),
    (["AAAAC", "AC"], Settings(k=3, unique=True), ("short records", "duplicates")),
    ([_MAIN] * 5 + [_BRANCH], Settings(k=4, clip_tips=True), ("tip clipping",)),
    ([random_dna(20, seed=1)], Settings(k=4, both_strands=True), ("both strands",)),
    (["ATGCATGCATGC"], Settings(k=3), ("cycle",)),
    (["ACGTNGGCA"], Settings(k=3), ("gc and n",)),
    (["LEMONLEMONADE"], Settings(k=3), ("letters",)),
)
_KINDS = {kind for *_, kinds in _CASES for kind in kinds}


def _analyses():
    return [
        analyze([Record(f"r{i}", s) for i, s in enumerate(reads, 1)], settings)
        for reads, settings, _ in _CASES
    ]


def _texts(a):
    texts = [a.headline, *a.notes, fasta_text(a)]
    texts += [header for table in TABLES for header in columns(table)]
    return texts + [genbank_text(a)] if a.dna else texts


def test_every_kind_of_note_is_given_by_some_input():
    # so that no note escapes the barred-word check by never appearing
    assert len(_KINDS) == 9
    for a, (_, _, kinds) in zip(_analyses(), _CASES, strict=True):
        assert len(a.notes) == len(kinds), a.notes


def test_no_output_names_a_count_as_coverage_or_depth():
    texts = [t for a in _analyses() for t in _texts(a)]
    assert [(i, barred_in(t)) for i, t in enumerate(texts) if barred_in(t)] == []


def test_gfa_and_csv_do_not_name_a_count_as_coverage_or_depth():
    texts = [t for a in _analyses() for t in (gfa_text(a), *(csv_text(a, n) for n in TABLES))]
    assert [(i, barred_in(t)) for i, t in enumerate(texts) if barred_in(t)] == []


@pytest.mark.parametrize(
    ("text", "found"),
    [
        (">ctg1 len=4 cov=10.0000x", ["cov="]),
        ("Removed 2 low-coverage edge(s)", ["coverage"]),
        ("Read DEPTH", ["depth"]),
        ("times seen; mean multiplicity", []),
    ],
)
def test_oracle_rejects_a_barred_word(text, found):
    assert barred_in(text) == found

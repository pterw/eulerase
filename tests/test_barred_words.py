"""No output calls an occurrence count "coverage" or "depth" (a check for barred words)."""

import pytest

from eulerase.analysis import TABLES, Settings, analyze, columns
from eulerase.io import fasta_text, genbank_text
from eulerase.sequence import Record

BARRED = ("coverage", "cov=", "depth")


def barred_in(text: str) -> list[str]:
    return [word for word in BARRED if word in text.lower()]


@pytest.mark.xfail(strict=True, reason="R3: coverage wording")
def test_no_output_names_a_count_as_coverage_or_depth():
    # GCA and CAT are seen twice, so a minimum multiplicity of 3 removes them and adds a note
    a = analyze([Record("read", "ATGCATGCATGC")], Settings(k=3, min_multiplicity=3))
    assert a.removed_kmers == 2 and a.notes
    texts = [fasta_text(a), genbank_text(a), *a.notes]
    texts += [header for table in TABLES for header in columns(table)]
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

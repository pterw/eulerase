"""Tests for the tables: one column definition, the number format and the row cap."""

import pytest
from genomes import random_dna

from eulerase.analysis import (
    MAX_TABLE_ROWS,
    NOT_AVAILABLE,
    TABLES,
    Analysis,
    Settings,
    analyze,
    columns,
    format_cell,
    limit_rows,
    rows,
    settings_text,
)
from eulerase.orfs import MIN_ORF_AA
from eulerase.sequence import Record

CONTIG_HEADERS = (
    "contig",
    "length",
    "edges",
    "k-mer count",
    "mean multiplicity",
    "min multiplicity",
    "gc",
    "orfs",
    "circular",
    "has n",
    "strand",
)


def _records(*sequences: str) -> list[Record]:
    return [Record(f"seq{i}", s) for i, s in enumerate(sequences, 1)]


def _yes_no(flag: bool) -> str:
    return "yes" if flag else "no"


def _fixed(x: float | None) -> str:
    return NOT_AVAILABLE if x is None else f"{x:.4f}"


def _expected(a: Analysis) -> dict[str, list[dict[str, str]]]:
    """The cells of every table, worked out from the analysis fields, keyed by header."""
    steps = [s for trail in a.trails for s in trail]
    m = a.metrics
    return {
        "Edges": [
            {
                "trail": str(s.trail),
                "step": str(s.step),
                "global step": str(s.global_step),
                "k-mer": s.kmer,
                "from": s.source,
                "to": s.target,
                "times seen": str(s.multiplicity),
            }
            for s in steps
        ],
        "Nodes": [
            {
                "node": n.node,
                "in": str(n.in_degree),
                "out": str(n.out_degree),
                "balance": str(n.balance),
                "role": n.role,
            }
            for n in a.nodes
        ],
        "Trails": [
            {"trail": str(i), "edges": str(len(t)), "sequence": seq}
            for i, (t, seq) in enumerate(zip(a.trails, a.spelled, strict=True), 1)
        ],
        "Contigs": [
            {
                "contig": str(c.id),
                "length": str(c.length),
                "edges": str(c.edges),
                "k-mer count": str(c.kmer_count),
                "mean multiplicity": _fixed(c.mean_multiplicity),
                "min multiplicity": str(c.min_multiplicity),
                "gc": _fixed(c.gc),
                "orfs": NOT_AVAILABLE if a.orfs is None else str(len(a.orfs[i])),
                "circular": _yes_no(c.circular),
                "has n": _yes_no(c.has_n),
                "strand": c.strand,
            }
            for i, c in enumerate(a.contigs)
        ],
        "Assembly": [
            {"measure": "contigs", "value": str(m.contigs)},
            {"measure": "total length", "value": str(m.total_length)},
            {"measure": "longest", "value": str(m.longest)},
            {"measure": "shortest", "value": str(m.shortest)},
            {"measure": "mean length", "value": _fixed(m.mean_length)},
            {"measure": "n50", "value": str(m.n50)},
            {"measure": "l50", "value": str(m.l50)},
            {"measure": "n90", "value": str(m.n90)},
            {"measure": "gc", "value": _fixed(m.gc)},
            {"measure": "mean multiplicity", "value": _fixed(m.mean_multiplicity)},
            {"measure": "distinct k-mers", "value": str(a.distinct_kmers)},
            {"measure": "total k-mers", "value": str(a.total_kmers)},
            {"measure": "removed k-mers", "value": str(a.removed_kmers)},
            {"measure": "clipped branches", "value": str(a.clipped_branches)},
            {"measure": "clipped k-mers", "value": str(a.clipped_kmers)},
        ],
    }


def _by_header(a: Analysis, table: str) -> list[dict[str, str]]:
    heads = columns(table)
    return [dict(zip(heads, row, strict=True)) for row in rows(a, table)]


def _matches(headers: tuple[str, ...], row: tuple[str, ...], expected: dict[str, str]) -> bool:
    return all(cell == expected[h] for h, cell in zip(headers, row, strict=True))


def test_oracle_rejects_a_shifted_column():
    expected = {"length": "6", "edges": "4"}
    row = ("6", "4")
    assert _matches(("length", "edges"), row, expected)
    assert not _matches(("length", "edges"), (row[1], row[0]), expected)


# The worked case: the read ATGCATGCATGC at k=3, by hand.
#   k-mers: ATG TGC GCA CAT ATG TGC GCA CAT ATG TGC, ten in all; ATG and TGC are seen 3 times,
#   GCA and CAT twice, so 4 distinct k-mers. They close a cycle AT-TG-GC-CA, one circular contig
#   ATGCAT: length 6 (k-1 plus 4 edges minus the 2 shared), edges 4, k-mer count 10,
#   mean multiplicity 10/4 = 2.5, minimum 2, GC 2/6 (G and C of ATGCAT) = 0.3333. It has no ORF
#   (no stop codon). One contig makes every length measure 6 and l50 1; mean length 6.
_WORKED_CONTIGS = [
    (
        "1",
        "6",
        "4",
        "10",
        "2.5000",
        "2",
        "0.3333",
        "0",
        "yes",
        "no",
        "as-assembled",
    )
]
_WORKED_ASSEMBLY = [
    ("contigs", "1"),
    ("total length", "6"),
    ("longest", "6"),
    ("shortest", "6"),
    ("mean length", "6.0000"),
    ("n50", "6"),
    ("l50", "1"),
    ("n90", "6"),
    ("gc", "0.3333"),
    ("mean multiplicity", "2.5000"),
    ("distinct k-mers", "4"),
    ("total k-mers", "10"),
    ("removed k-mers", "0"),
    ("clipped branches", "0"),
    ("clipped k-mers", "0"),
]


@pytest.mark.parametrize(
    ("reads", "settings", "worked", "differs"),
    [
        (["ATGCATGCATGC"], Settings(k=3), True, None),
        (["ATGCATGCATGC", "GATTACAGATTACA", "CCCGGGTTTAAAC"], Settings(k=3), False, None),
        # Contigs of 36 and 30 letters: n50 (first contig) differs from n90 (both needed).
        (
            ["ATGGCTGCCGCAGCGCTGCTACTCCTTGGTGGCTAA", "ACGTTGCAAGGCTTAGCCATCCGTACCAGT"],
            Settings(k=12),
            False,
            "n50",
        ),
        # One dead-end branch of three k-mers seen once, beside a path seen five times.
        (
            ["AGGAGTTAAATCGATGTCTCCTTCTGGCTT"] * 5 + ["AGGAGTTAAATCGATGCCA"],
            Settings(k=4, clip_tips=True),
            False,
            "clipped",
        ),
    ],
)
def test_every_cell_equals_a_value_computed_from_the_analysis(reads, settings, worked, differs):
    a = analyze(_records(*reads), settings)
    if differs == "n50":
        assert [c.length for c in a.contigs] == [36, 30]
        assert a.metrics.n50 != a.metrics.n90
    if differs == "clipped":
        assert a.clipped_branches != a.clipped_kmers
        assert a.clipped_branches > 0
        assert a.clipped_kmers > 0
    expected = _expected(a)
    for table in TABLES:
        got = rows(a, table)
        want = expected[table]
        extra = 2 if table == "Assembly" else 0  # the settings and orf rule rows
        assert len(got) == len(want) + extra
        heads = columns(table)
        assert all(_matches(heads, row, w) for row, w in zip(got[: len(want)], want, strict=True))
        assert set(want[0]) == set(heads)
    assert columns("Contigs") == CONTIG_HEADERS
    if worked:
        assert rows(a, "Contigs") == _WORKED_CONTIGS
        assert rows(a, "Assembly")[: len(_WORKED_ASSEMBLY)] == _WORKED_ASSEMBLY


@pytest.mark.parametrize(
    ("value", "text"),
    [
        (7, "7"),
        (1234567, "1234567"),
        (0.5, "0.5000"),
        (2 / 3, "0.6667"),
        (True, "yes"),
        (False, "no"),
        (None, "not available"),
        ("x", "x"),
    ],
)
def test_format_cell_cases(value, text):
    assert format_cell(value) == text


def test_gc_none_shows_not_available():
    a = analyze(_records("NNNNNN"), Settings(k=3))
    assert a.contigs[0].gc is None
    assert _by_header(a, "Contigs")[0]["gc"] == "not available"


def test_orfs_column_is_not_available_for_non_dna():
    a = analyze(_records("MKVLAAGIV"), Settings(k=3))
    assert a.orfs is None
    assert {r["orfs"] for r in _by_header(a, "Contigs")} == {"not available"}


def test_limit_rows_says_how_many_of_how_many():
    a = analyze(_records(random_dna(2502, seed=7)), Settings(k=3))
    all_rows = rows(a, "Edges")
    assert len(all_rows) == 2500
    cut, note = limit_rows(all_rows)
    assert cut == all_rows[:2000]
    assert note == "Showing 2000 of 2500 rows."
    assert limit_rows(all_rows[:MAX_TABLE_ROWS]) == (all_rows[:2000], None)


def test_assembly_table_lists_the_settings_and_the_removals():
    settings = Settings(k=3, min_multiplicity=2)
    a = analyze(_records("ATGCATGCATGC", "GATTACA"), settings)
    table = dict(rows(a, "Assembly"))
    assert a.removed_kmers > 0
    assert table["settings"] == settings_text(settings)
    assert table["removed k-mers"] == str(a.removed_kmers)
    rule = table["orf rule"]
    assert "ATG" in rule
    assert "standard" in rule
    assert str(MIN_ORF_AA) in rule
    protein = analyze(_records("MKVLAAGIV"), Settings(k=3))
    assert dict(rows(protein, "Assembly"))["orf rule"] == "not available"

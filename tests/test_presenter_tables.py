"""The presenter's table views: every cell sits under the header it belongs to."""

import random

import pytest
from genomes import random_dna

from eulerase import analysis
from eulerase.analysis import NOT_AVAILABLE, Settings, analyze
from eulerase.gui.presenter import records_to_text, result_view, table_views
from eulerase.sequence import Record, parse_records

READS = [
    Record("a", random_dna(30, 1) + "ACGTACG" + random_dna(10, 2) + "ACGTACG"),
    Record("b", random_dna(25, 3) + "ACGTACG" + random_dna(8, 4)),
]


def dna_analysis(**settings):
    return analyze(READS, Settings(k=7, **settings))


def plain(value):
    """The formatting the spec states, written out again here."""
    if value is None:
        return NOT_AVAILABLE
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def contig_expected(contig, orf_count):
    return {
        "contig": contig.id,
        "length": contig.length,
        "edges": contig.edges,
        "k-mer count": contig.kmer_count,
        "mean multiplicity": contig.mean_multiplicity,
        "min multiplicity": contig.min_multiplicity,
        "gc": contig.gc,
        "orfs": orf_count,
        "cycle": contig.circular,
        "has n": contig.has_n,
        "strand": contig.strand,
    }


def edge_expected(step):
    return {
        "trail": step.trail,
        "step": step.step,
        "global step": step.global_step,
        "k-mer": step.kmer,
        "from": step.source,
        "to": step.target,
        "times seen": step.multiplicity,
    }


def by_header(table):
    return [dict(zip(table.columns, row, strict=True)) for row in table.rows]


def expected_cells(a):
    contigs = [contig_expected(c, len(o)) for c, o in zip(a.contigs, a.orfs, strict=True)]
    edges = [edge_expected(s) for trail in a.trails for s in trail]
    return {
        "Contigs": [{h: plain(v) for h, v in row.items()} for row in contigs],
        "Edges": [{h: plain(v) for h, v in row.items()} for row in edges],
    }


def cells_match(a, views):
    shown = {t.name: by_header(t) for t in views}
    return all(shown[name] == rows for name, rows in expected_cells(a).items())


def test_every_cell_under_a_header_equals_a_value_computed_from_the_contigs():
    a = dna_analysis()
    assert a.contigs and a.orfs is not None
    views = table_views(a)
    assert [t.name for t in views] == list(analysis.TABLES)
    assert cells_match(a, views)


def test_oracle_rejects_a_shifted_column():
    a = dna_analysis()
    views = list(table_views(a))
    contigs = next(t for t in views if t.name == "Contigs")
    row = list(contigs.rows[0])
    assert row[1] != row[2]  # length and edges differ, so a swap changes the row
    row[1], row[2] = row[2], row[1]
    shifted = type(contigs)(contigs.name, contigs.columns, (tuple(row), *contigs.rows[1:]), None)
    views[views.index(contigs)] = shifted
    assert not cells_match(a, views)


@pytest.mark.parametrize("table", analysis.TABLES)
@pytest.mark.parametrize(
    "records, settings",
    [
        (READS, Settings(k=7)),
        (parse_records("LEMONLEMONADE"), Settings(k=3)),
        (READS, Settings(k=7, min_multiplicity=99)),
    ],
    ids=["dna", "text", "contigs emptied by min multiplicity"],
)
def test_every_row_has_as_many_cells_as_headers_in_every_table(table, records, settings):
    view = next(t for t in table_views(analyze(records, settings)) if t.name == table)
    assert len(view.columns) > 0
    assert all(len(row) == len(view.columns) for row in view.rows)


def test_a_cut_table_says_how_many_rows_it_shows():
    a = analyze([Record("long", random_dna(analysis.MAX_TABLE_ROWS + 60, 5))], Settings(k=12))
    total = len(analysis.rows(a, "Edges"))
    assert total > analysis.MAX_TABLE_ROWS
    views = {t.name: t for t in table_views(a)}
    assert len(views["Edges"].rows) == analysis.MAX_TABLE_ROWS
    assert views["Edges"].note is not None
    assert str(total) in views["Edges"].note
    short = table_views(dna_analysis())
    assert all(t.note is None for t in short if t.name != "Contigs")


@pytest.mark.parametrize("cut", [False, True], ids=["whole table", "cut table"])
def test_the_contigs_table_states_the_orf_rule_after_any_cut_note(cut, monkeypatch):
    a = dna_analysis()
    total = len(a.contigs)
    assert total > 0
    if cut:
        monkeypatch.setattr(analysis, "MAX_TABLE_ROWS", total - 1)
    rule_line = f"orf rule: {dict(analysis.rows(a, 'Assembly'))['orf rule']}"
    view = result_view(a, None)
    notes = {t.name: t.note for t in view.tables}
    cut_line = f"Showing {total - 1} of {total} rows.\n"
    assert notes["Contigs"] == (cut_line if cut else "") + rule_line
    assert rule_line in view.details
    assert all("orf rule" not in (note or "") for name, note in notes.items() if name != "Contigs")


def test_settings_and_every_removal_note_reach_the_result():
    a = analyze(READS, Settings(k=7, min_multiplicity=2, clip_tips=True))
    assert a.notes
    view = result_view(a, "layout")
    assert view.details[0] == f"settings: {analysis.settings_text(a.settings)}"
    assert view.details[1].startswith("orf rule: ")
    assert all(note in view.details for note in a.notes)
    assert a.headline in view.status
    assert analysis.settings_text(a.settings) in view.status
    assert view.steps == sum(len(trail) for trail in a.trails)
    assert view.layout == "layout"


def test_the_assembly_table_holds_only_numeric_rows():
    view = next(t for t in table_views(dna_analysis()) if t.name == "Assembly")
    measures = [row[0] for row in view.rows]
    assert measures[0] == "contigs" and measures[-1] == "clipped k-mers"
    assert "settings" not in measures and "orf rule" not in measures


@pytest.mark.parametrize("seed", range(4))
def test_records_to_text_round_trips_random_records(seed):
    lengths = random.Random(seed).sample(range(1, 60), 3)
    records = [Record(f"r{i}", random_dna(n, seed * 10 + i)) for i, n in enumerate(lengths, 1)]
    assert parse_records(records_to_text(records)) == records

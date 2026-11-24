"""Loading a file and exporting: what is asked, what is written, and what stays untouched."""

from dataclasses import replace

import pytest
from genomes import random_dna
from gui_fakes import SHORT_READ, FakeScheduler, FakeView, FakeWorker

from eulerase import io, kmers
from eulerase.analysis import MAX_TABLE_ROWS, TABLES, Settings, analyze
from eulerase.gui.presenter import BUSY, EMPTY, READY, FormValues, Presenter
from eulerase.sequence import Record, parse_records

OVERLAPPING = ">r1\nACGTGCAATG\n>r2\nCAATGCCGTA\n>r3\nCCGTAACGTG"
FIVE_FILES = ("x_edges.csv", "x_nodes.csv", "x_trails.csv", "x_contigs.csv", "x_assembly.csv")
CONTIG_EXPORTS = ("fasta", "genbank", "gfa")
WRITERS = {"fasta": io.fasta_text, "genbank": io.genbank_text, "gfa": io.gfa_text}


def form(text, k_text="5", min_multiplicity_text=""):
    return FormValues(text, "sequences", k_text, False, False, min_multiplicity_text, False)


def rig(start=SHORT_READ):
    view = FakeView(start)
    worker = FakeWorker(view.calls)
    return view, worker, Presenter(view, worker, FakeScheduler())


def built(start):
    view, worker, presenter = rig(start)
    presenter.build()
    worker.finish(0, presenter)
    assert presenter.phase == READY
    return view, presenter


def settings_of(start):
    return Settings(
        k=int(start.k_text), min_multiplicity=int(start.min_multiplicity_text or 0) or None
    )


def write_fasta(path, records):
    path.write_text("".join(f">{r.name}\n{r.sequence}\n" for r in records))


def dialog_lists(calls, paths, existing):
    """Whether the log holds one export dialog listing exactly these paths and existing ones."""
    asked = [call for call in calls if call[0] == "confirm_export"]
    return asked == [("confirm_export", tuple(paths), tuple(existing))]


def asks_nothing_and_writes_nothing(view, tmp_path, action):
    calls, files = list(view.calls), sorted(tmp_path.iterdir())
    action()
    assert view.calls == calls
    assert sorted(tmp_path.iterdir()) == files


@pytest.mark.parametrize("kind", CONTIG_EXPORTS)
def test_each_contig_export_writes_what_the_core_writer_returns(kind, tmp_path):
    start = form(OVERLAPPING)
    view, presenter = built(start)
    expected = WRITERS[kind](analyze(parse_records(OVERLAPPING), settings_of(start)))
    assert expected  # the writer has something to say, so an empty file would differ
    view.save_path = str(tmp_path / "out.txt")

    presenter.export(kind)

    data = (tmp_path / "out.txt").read_bytes()
    assert data == expected.encode("utf-8")
    assert b"\r" not in data


def test_the_figure_export_hands_the_analysis_and_the_chosen_path_to_the_window(tmp_path):
    view, presenter = built(form(OVERLAPPING))
    view.save_path = str(tmp_path / "pic.png")

    presenter.export("figure")

    assert ("save_figure", presenter.analysis, view.save_path) in view.calls
    assert list(tmp_path.iterdir()) == []  # the window writes the picture, not the presenter


@pytest.mark.parametrize(
    "typed, present",
    [("x.csv", ("x_edges.csv", "x_nodes.csv")), ("x.csv", ()), ("x.CSV", ())],
    ids=["two files exist", "none exist", "upper-case suffix"],
)
def test_the_csv_dialog_lists_every_file_it_will_write_and_asks_before_replacing_any(
    typed, present, tmp_path
):
    view, presenter = built(form(OVERLAPPING))
    for name in present:
        (tmp_path / name).write_text("old")
    view.save_path = str(tmp_path / typed)
    expected = [str(tmp_path / name) for name in FIVE_FILES]

    presenter.export("csv")

    assert dialog_lists(view.calls, expected, [str(tmp_path / n) for n in present])


def test_oracle_rejects_a_dialog_that_omits_a_file(tmp_path):
    paths = [str(tmp_path / name) for name in FIVE_FILES]
    complete, omitting = FakeView(), FakeView()
    complete.confirm_export(paths, [])
    omitting.confirm_export(paths[:4], [])
    assert dialog_lists(complete.calls, paths, [])  # the check accepts the right dialog
    assert not dialog_lists(omitting.calls, paths, [])


def test_declining_the_csv_dialog_writes_nothing_and_leaves_existing_files_unchanged(tmp_path):
    view, presenter = built(form(OVERLAPPING))
    (tmp_path / "x_edges.csv").write_text("old")
    view.save_path, view.export_answer = str(tmp_path / "x.csv"), False

    presenter.export("csv")

    assert sorted(p.name for p in tmp_path.iterdir()) == ["x_edges.csv"]
    assert (tmp_path / "x_edges.csv").read_text() == "old"


def test_accepting_the_csv_dialog_writes_all_five_tables_with_every_row(tmp_path):
    start = form(random_dna(2200, 5), k_text="12")
    view, presenter = built(start)
    expected = analyze(parse_records(start.text), settings_of(start))
    edges = io.csv_text(expected, "Edges")
    assert edges.count("\n") - 1 > MAX_TABLE_ROWS  # more rows than the table view shows
    view.save_path = str(tmp_path / "x.csv")

    presenter.export("csv")

    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(FIVE_FILES)
    for table in TABLES:
        written = (tmp_path / f"x_{table.lower()}.csv").read_bytes()
        assert written == io.csv_text(expected, table).encode("utf-8")


def test_export_is_unavailable_until_a_build_succeeds(tmp_path):
    view, worker, presenter = rig()
    view.save_path = str(tmp_path / "x.csv")
    kinds = ("figure", "csv", *CONTIG_EXPORTS)

    def export_all():
        for kind in kinds:
            presenter.export(kind)

    asks_nothing_and_writes_nothing(view, tmp_path, export_all)  # before any build
    presenter.build()
    asks_nothing_and_writes_nothing(view, tmp_path, export_all)  # while busy
    worker.finish(0, presenter)
    view.form = replace(SHORT_READ, k_text="x")
    presenter.build()
    assert presenter.phase == EMPTY
    asks_nothing_and_writes_nothing(view, tmp_path, export_all)  # after a failed build


def test_contig_exports_are_unavailable_when_no_contig_was_assembled(tmp_path):
    start = form("ACGTACGGATTGCCA", k_text="4", min_multiplicity_text="2")
    view, presenter = built(start)
    assert presenter.analysis.contigs == ()  # the minimum removed every k-mer from assembly
    view.save_path = str(tmp_path / "out.txt")

    def export_contigs():
        for kind in CONTIG_EXPORTS:
            presenter.export(kind)

    asks_nothing_and_writes_nothing(view, tmp_path, export_contigs)


@pytest.mark.parametrize(
    "kind, saved, blocker, failed",
    [
        ("csv", "missing/x.csv", None, "x_edges.csv"),
        ("fasta", "missing/out.fasta", None, "out.fasta"),
        ("csv", "x.csv", "x_nodes.csv", "x_nodes.csv"),  # the second file of the five
    ],
    ids=["csv, no folder", "fasta, no folder", "csv, second file blocked"],
)
def test_a_write_failure_is_one_sentence_naming_the_file_that_failed(
    kind, saved, blocker, failed, tmp_path
):
    view, presenter = built(form(OVERLAPPING))
    view.save_path = str(tmp_path / saved)
    if blocker:
        (tmp_path / blocker).mkdir()  # a directory cannot be opened for writing

    presenter.export(kind)

    assert view.status_is_error
    assert view.status and "\n" not in view.status and "Traceback" not in view.status
    assert failed in view.status
    assert presenter.phase == READY  # the result stays; only the write failed
    if kind == "csv":
        written = sum((tmp_path / name).is_file() for name in FIVE_FILES)
        assert written != len(FIVE_FILES)
        # how many of the set were written before it (the folder's own name may hold digits)
        assert str(written) in view.status.replace(str(tmp_path), "")


@pytest.mark.parametrize(
    "ask, action",
    [
        ("ask_open_path", lambda p: p.open_file()),
        ("ask_save_path", lambda p: p.export("figure")),
        ("ask_save_path", lambda p: p.export("fasta")),
        ("ask_save_path", lambda p: p.export("csv")),
    ],
    ids=["open", "figure", "fasta", "csv"],
)
def test_cancelling_a_dialog_changes_nothing(ask, action, tmp_path):
    view, worker, presenter = rig(form(OVERLAPPING))
    presenter.build()
    worker.finish(0, presenter)
    shown, files, mark = view.result, sorted(tmp_path.iterdir()), len(view.calls)
    view.open_path = view.save_path = None

    action(presenter)

    assert [call[0] for call in view.calls[mark:]] == [ask]  # asked, and nothing after
    assert (view.result, presenter.phase) == (shown, READY)
    assert len(worker.jobs) == 1
    assert sorted(tmp_path.iterdir()) == files


def test_loading_a_file_shows_every_record_in_the_input_and_builds(tmp_path):
    records = [Record("a", random_dna(40_000, 1)), Record("b", random_dna(30_000, 2))]
    assert len(records[0].sequence) != len(records[1].sequence)
    assert sum(len(r.sequence) for r in records) > 60_000
    path = tmp_path / "reads.fasta"
    write_fasta(path, records)
    view, worker, presenter = rig(replace(SHORT_READ, input_mode="kmers", k_text="4"))
    presenter.mode_changed("kmers")
    assert view.k_enabled is False

    view.open_path = str(path)
    presenter.open_file()

    assert [(r.name, r.sequence) for r in parse_records(view.form.text)] == [
        (r.name, r.sequence) for r in records
    ]
    assert view.form.input_mode == "sequences"
    assert view.form.k_text == "4"  # the other controls are unchanged
    assert view.k_enabled is True
    assert len(worker.jobs) == 1


def test_a_file_past_the_k_mer_limit_is_loaded_whole_and_the_build_refuses_it(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(kmers, "MAX_KMERS", 1000)
    records = [Record("a", random_dna(3000, 3))]
    path = tmp_path / "long.fasta"
    write_fasta(path, records)
    view, worker, presenter = rig(replace(SHORT_READ, k_text="4"))
    view.open_path = str(path)

    presenter.open_file()

    assert [(r.name, r.sequence) for r in parse_records(view.form.text)] == [
        ("a", records[0].sequence)
    ]
    assert presenter.phase == BUSY and len(worker.jobs) == 1
    worker.finish(0, presenter)
    assert presenter.phase == EMPTY
    assert view.status_is_error
    assert f"{kmers.MAX_KMERS:,}" in view.status and "\n" not in view.status


@pytest.mark.parametrize("kind", ["missing", "directory", "empty"])
def test_a_file_that_cannot_be_loaded_leaves_the_screen_alone(kind, tmp_path):
    view, worker, presenter = rig(form(OVERLAPPING))
    presenter.build()
    worker.finish(0, presenter)
    shown, start, mark = view.result, view.form, len(view.calls)
    path = tmp_path / "input"
    if kind == "directory":
        path.mkdir()
    elif kind == "empty":
        path.write_bytes(b"")

    presenter.load_file(str(path))

    assert view.status_is_error
    assert view.status and "\n" not in view.status
    assert (view.form, view.result, presenter.phase) == (start, shown, READY)
    assert len(worker.jobs) == 1  # no new build
    assert view.calls[mark:] == [("show_status", view.status, True)]  # nothing else was touched

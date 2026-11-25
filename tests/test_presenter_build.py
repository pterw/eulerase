"""The build lifecycle: a build runs off the window's thread and ends busy, stopped or failed."""

from dataclasses import replace

import pytest
from gui_fakes import SHORT_READ, FakeScheduler, FakeView, FakeWorker, protocol_gaps

from eulerase.analysis import Settings, analyze
from eulerase.gui.presenter import (
    BUSY,
    EMPTY,
    EXPORT_KINDS,
    READY,
    Presenter,
    Scheduler,
    View,
    Worker,
)
from eulerase.sequence import parse_records


def rig(form=SHORT_READ, layout=None):
    view = FakeView(form)
    worker = FakeWorker(view.calls)
    return view, worker, Presenter(view, worker, FakeScheduler(), layout=layout)


def results_shown_since(calls, mark):
    """The history check: every result put on screen after position `mark`."""
    return [call for call in calls[mark:] if call[0] == "show_result"]


def test_a_build_runs_in_the_worker_with_a_busy_state_and_shows_the_analysis():
    view, worker, presenter = rig()
    presenter.build()
    assert len(worker.jobs) == 1
    assert (view.busy, view.result, presenter.analysis) == (True, None, None)
    assert presenter.phase == BUSY

    worker.finish(0, presenter)
    expected = analyze(parse_records(SHORT_READ.text), Settings(k=3))
    assert view.result.analysis == expected
    assert presenter.analysis == expected
    assert (view.busy, presenter.phase) == (False, READY)
    assert presenter.step == presenter.total == sum(len(t) for t in expected.trails)
    assert presenter.total > 0
    assert ("show_step", presenter.total, presenter.total) in view.calls
    assert ("set_playing", False) in view.calls
    assert ("show_status", view.result.status, False) in view.calls

    first_status = view.status
    assert first_status  # so a difference below cannot come from an empty line
    presenter.build()
    assert presenter.phase == BUSY
    assert view.status and view.status != first_status


@pytest.mark.parametrize("phase", [EMPTY, READY])
def test_stop_does_nothing_unless_a_build_is_running(phase):
    view, worker, presenter = rig()
    if phase == READY:
        presenter.build()
        worker.finish(0, presenter)
    calls, result, exports = list(view.calls), view.result, dict(view.exports)

    presenter.stop()

    assert presenter.phase == phase
    assert (view.calls, view.result, view.exports) == (calls, result, exports)


@pytest.mark.parametrize(
    "form, abandon",
    [(SHORT_READ, Presenter.stop), (replace(SHORT_READ, k_text="abc"), Presenter.build)],
    ids=["stop", "a build that fails"],
)
def test_an_abandoned_build_is_discarded_when_its_result_arrives(form, abandon):
    view, worker, presenter = rig()
    presenter.build()
    assert worker.run_without_delivering(0).analysis is not None  # it would have succeeded
    mark = len(view.calls)

    view.form = form
    abandon(presenter)
    worker.finish(0, presenter)

    assert results_shown_since(view.calls, mark) == []
    assert (presenter.phase, view.busy) == (EMPTY, False)
    assert view.exports == dict.fromkeys(EXPORT_KINDS, False)


@pytest.mark.parametrize("order", [(0, 1), (1, 0)], ids=["old first", "old last"])
def test_a_new_build_abandons_the_one_in_progress(order):
    view, worker, presenter = rig()
    presenter.build()
    view.form = replace(SHORT_READ, k_text="4")
    presenter.build()
    assert len(worker.jobs) == 2

    for i in order:
        worker.finish(i, presenter)

    shown = results_shown_since(view.calls, 0)
    assert [call[1].analysis.k for call in shown] == [4]
    assert (presenter.phase, presenter.analysis.k) == (READY, 4)


def good_then(**changes):
    return replace(SHORT_READ, **changes)


FAILURES = [
    pytest.param(good_then(text=""), None, id="empty input"),
    pytest.param(good_then(text="LEMONLEMON", both_strands=True), None, id="both strands on text"),
    pytest.param(good_then(k_text="99"), None, id="k above every record"),
    pytest.param(good_then(k_text="abc"), None, id="k not a number"),
    pytest.param(good_then(k_text="4.0"), None, id="k with a decimal point"),
    pytest.param(good_then(k_text="1"), None, id="k below 2"),
    pytest.param(good_then(min_multiplicity_text="0"), None, id="minimum below 1"),
    pytest.param(good_then(min_multiplicity_text="x"), None, id="minimum not a number"),
    pytest.param(SHORT_READ, RuntimeError("boom"), id="a bug in the job"),
]
REFUSED_BEFORE_ANY_JOB = {"k below 2"}  # an invalid k must never reach the worker


@pytest.mark.parametrize("bad_form, bug", FAILURES)
def test_a_failed_build_clears_the_result_and_export_is_unavailable(request, bad_form, bug):
    armed = []

    def layout(a, aspect):
        if armed:
            raise bug

    view, worker, presenter = rig(layout=layout)
    presenter.build()
    worker.finish(0, presenter)
    assert presenter.phase == READY and any(view.exports.values())  # something to clear

    view.form = bad_form
    armed.extend([bug] if bug else [])
    mark = len(view.calls)
    presenter.build()
    if request.node.callspec.id in REFUSED_BEFORE_ANY_JOB:
        assert len(worker.jobs) == 1  # only the good build's job
    if len(worker.jobs) > 1:  # this failure is found inside the job
        worker.finish(1, presenter)

    assert presenter.phase == EMPTY
    assert presenter.analysis is None
    assert ("clear_result",) in view.calls[mark:]
    assert view.exports == dict.fromkeys(EXPORT_KINDS, False)
    assert view.busy is False
    assert view.status_is_error
    assert view.status and "\n" not in view.status and "Traceback" not in view.status
    assert type(bug).__name__ in view.status if bug else True


@pytest.mark.parametrize(
    "minimum", ["", "   ", "9"], ids=["contigs", "only spaces means none", "none left"]
)
def test_contig_exports_are_enabled_only_when_the_analysis_has_contigs(minimum):
    view, worker, presenter = rig(good_then(min_multiplicity_text=minimum))
    presenter.build()
    worker.finish(0, presenter)
    settings = Settings(k=3, min_multiplicity=int(minimum) if minimum.strip() else None)
    has_contigs = bool(analyze(parse_records(SHORT_READ.text), settings).contigs)
    assert has_contigs == (minimum != "9")  # the cases differ in the thing tested
    assert view.exports == {
        "figure": True,
        "csv": True,
        "fasta": has_contigs,
        "genbank": has_contigs,
        "gfa": has_contigs,
    }


def test_a_control_that_cannot_be_edited_is_never_validated():
    view, worker, presenter = rig(
        replace(SHORT_READ, text="ACG CGT GTA TAC", input_mode="kmers", k_text="abc")
    )
    presenter.build()
    worker.finish(0, presenter)
    assert presenter.phase == READY


def test_k_typed_with_spaces_builds():
    view, worker, presenter = rig(good_then(k_text=" 5 "))
    presenter.build()
    worker.finish(0, presenter)
    assert presenter.analysis.k == 5
    assert presenter.analysis.k != Settings().k  # not the default by accident


def test_the_layout_receives_the_aspect_read_with_the_form():
    asked = []
    view, worker, presenter = rig(layout=lambda a, aspect: asked.append((a, aspect)))
    first, later = 1.7, 0.4
    view.aspect = first
    presenter.build()
    view.aspect = later  # the window is resized before the job runs
    worker.finish(0, presenter)

    assert asked == [(analyze(parse_records(SHORT_READ.text), Settings(k=3)), first)]
    reads = [call[0] for call in view.calls]
    assert reads.index("read_form") < reads.index("canvas_aspect") < reads.index("worker_start")


@pytest.mark.parametrize("mode", ["sequences", "kmers"])
def test_mode_change_enables_k_only_for_sequences(mode):
    view, _, presenter = rig()
    assert view.k_enabled is None
    presenter.mode_changed(mode)
    assert view.k_enabled is (mode == "sequences")


@pytest.mark.parametrize(
    "fake, protocol",
    [(FakeView, View), (FakeWorker, Worker), (FakeScheduler, Scheduler)],
    ids=lambda x: x.__name__,
)
def test_the_fakes_implement_the_protocols(fake, protocol):
    assert protocol_gaps(fake, protocol) == []


def test_oracle_rejects_a_class_missing_a_protocol_method():
    class RenamedParameter(FakeView):
        def show_status(self, text, error=False): ...

    kept = {n: v for n, v in vars(FakeView).items() if n != "show_step"}
    missing = type("MissingShowStep", (), kept)
    assert [g.split(":")[0] for g in protocol_gaps(missing, View)] == ["show_step"]
    assert [g.split(":")[0] for g in protocol_gaps(RenamedParameter, View)] == ["show_status"]


def test_oracle_rejects_a_result_shown_after_stop():
    forged = [("show_busy", True), ("show_status", "stopped", False), ("show_result", object())]
    assert len(results_shown_since(forged, 2)) == 1
    assert results_shown_since(forged, 3) == []

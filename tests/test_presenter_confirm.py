"""The confirmation before a large build: asked once, in the stated range; a no changes nothing.
A build the question covers runs in a process, any other in the thread."""

import pickle
from dataclasses import replace

import pytest
from gui_fakes import SHORT_READ, FakeScheduler, FakeView, FakeWorker

from eulerase import analysis, render, sequence
from eulerase.gui.presenter import (
    BUSY,
    EMPTY,
    READY,
    FormValues,
    Presenter,
    run_build,
)

THRESHOLD = 50
K = 5


def big_form(length: int, both_strands: bool = False) -> FormValues:
    """One sequence of `length` letters: it holds `length - 4` 5-mers on one strand."""
    text = ("ACGT" * (length // 4 + 1))[:length]
    return FormValues(text, "sequences", str(K), both_strands, False, "", False)


@pytest.fixture(autouse=True)
def low_threshold(monkeypatch):
    monkeypatch.setattr("eulerase.gui.presenter.CONFIRM_ABOVE_KMERS", THRESHOLD)


def rig():
    view = FakeView(SHORT_READ)
    scheduler = FakeScheduler()
    worker = FakeWorker(view.calls)
    return view, worker, scheduler, Presenter(view, worker, scheduler, layout=render.compute_layout)


def confirms(view):
    return [call for call in view.calls if call[0] == "confirm"]


@pytest.mark.parametrize(
    ("length", "both_strands", "kmer_count", "asks", "limit"),
    [
        (54, False, 50, False, None),  # equal to the threshold: no question
        (55, False, 51, True, None),  # one above it: asks
        (29, True, 50, False, None),  # 25 k-mers, counted twice, equal to the threshold
        (30, True, 52, True, None),  # 26 k-mers, counted twice, above it
        (1004, False, 1000, True, None),  # the question names the estimate, thousands separated
        (55, False, 51, True, 51),  # equal to the core limit: still asked, then built
    ],
)
def test_a_build_asks_only_above_the_threshold_counting_both_strands_twice(
    length, both_strands, kmer_count, asks, limit, monkeypatch
):
    if limit:
        monkeypatch.setattr("eulerase.kmers.MAX_KMERS", limit)
        assert limit == kmer_count  # exactly at the limit, not below it
    view, worker, _, presenter = rig()
    view.form = big_form(length, both_strands)
    # By hand: (length - K + 1) k-mers on a strand, doubled with both strands.
    assert kmer_count == (length - K + 1) * (2 if both_strands else 1)
    assert (kmer_count > THRESHOLD) == asks  # the cases sit on each side of the line

    presenter.build()

    assert len(confirms(view)) == int(asks)
    if asks:
        assert f"{kmer_count:,}" in confirms(view)[0][1]
    # a build the question covers goes to a process, every other to the thread
    assert (len(worker.jobs), len(worker.processes), presenter.phase) == (
        (0, 1, BUSY) if asks else (1, 0, BUSY)
    )


def test_a_large_build_hands_its_process_a_picklable_run_build_giving_the_threads_outcome():
    view, worker, _, presenter = rig()
    view.form = big_form(55)
    presenter.build()
    ((token, function, args),) = worker.processes
    pickle.dumps(args)  # a closure, or the default lambda layout, would not pickle
    assert function is run_build

    outcome = function(token, *args)

    records = sequence.parse_records(view.form.text)
    expected = analysis.analyze(records, analysis.Settings(k=K))
    assert (outcome.analysis, outcome.error) == (expected, None)
    assert outcome.layout == render.compute_layout(outcome.analysis, view.aspect)
    presenter.on_build_outcome(outcome)
    assert presenter.phase == READY, "the token is the one the build took"


@pytest.mark.parametrize("state", ["walk playing", "build running"])
def test_declining_leaves_the_screen_the_build_and_the_walk_as_they_were(state):
    view, worker, scheduler, presenter = rig()
    presenter.build()
    if state == "walk playing":
        worker.finish(0, presenter)
        presenter.set_step(0)
        presenter.play()
    view.form = big_form(55)
    view.confirm_answer = False
    before = (presenter.phase, presenter.analysis, view.result, dict(scheduler.pending))
    assert len(scheduler.pending) == (1 if state == "walk playing" else 0)
    mark, jobs, processes = len(view.calls), len(worker.jobs), len(worker.processes)

    presenter.build()

    names = [call[0] for call in view.calls[mark:]]
    assert names == ["read_form", "canvas_aspect", "confirm", "show_status"]
    assert view.status_is_error is False
    assert (presenter.phase, presenter.analysis, view.result, scheduler.pending) == before
    assert (len(worker.jobs), len(worker.processes)) == (jobs, processes)
    assert worker.ended == 0, "a declined question abandons nothing"
    assert presenter.playing == (state == "walk playing")
    if state == "build running":  # not abandoned: its outcome still arrives
        worker.finish(0, presenter)
        assert (presenter.phase, view.result is not None) == (READY, True)


def test_an_input_above_the_core_limit_is_refused_not_confirmed(monkeypatch):
    monkeypatch.setattr("eulerase.kmers.MAX_KMERS", 40)
    view, worker, _, presenter = rig()
    view.form = big_form(30, True)  # 52 k-mers: above the threshold and the limit

    presenter.build()
    assert (len(worker.jobs), len(worker.processes)) == (1, 0)  # the thread way, as before
    worker.finish(0, presenter)

    assert confirms(view) == []
    assert (presenter.phase, view.status_is_error) == (EMPTY, True)
    assert "40" in view.status


def running_large_build():
    """A first large build in flight: it started a process and ended none."""
    view, worker, scheduler, presenter = rig()
    view.form = big_form(55)
    presenter.build()
    assert (len(worker.processes), worker.ended, presenter.phase) == (1, 0, BUSY)
    return view, worker, scheduler, presenter


def test_stop_ends_the_process_and_the_abandoned_outcome_changes_nothing():
    view, worker, _, presenter = running_large_build()

    presenter.stop()

    assert (worker.ended, presenter.phase) == (1, EMPTY)
    mark = len(view.calls)
    token, function, args = worker.processes[0]
    presenter.on_build_outcome(function(token, *args))  # as if it had finished just before
    assert (presenter.phase, view.result, view.calls[mark:]) == (EMPTY, None, [])


@pytest.mark.parametrize("before", ["no build yet", "a finished large build"])
def test_stop_with_no_build_running_ends_nothing(before):
    view, worker, _, presenter = rig()
    if before == "a finished large build":
        view.form = big_form(55)
        presenter.build()
        token, function, args = worker.processes[0]
        presenter.on_build_outcome(function(token, *args))
        assert presenter.phase == READY

    presenter.stop()

    assert worker.ended == 0


@pytest.mark.parametrize(
    ("first", "new_length", "ends"),
    [
        ("running", 55, 1),  # a new large build
        ("running", 20, 1),  # a new small build: the running process ends all the same
        ("finished", 55, 0),  # nothing is running, so nothing is ended
    ],
)
def test_a_new_build_ends_a_running_process_before_it_starts_its_own(first, new_length, ends):
    view, worker, _, presenter = running_large_build()
    if first == "finished":
        token, function, args = worker.processes[0]
        presenter.on_build_outcome(function(token, *args))
        assert presenter.phase == READY
    view.form = big_form(new_length)
    mark = len(view.calls)

    presenter.build()

    assert worker.ended == ends
    new = "worker_start_process" if new_length == 55 else "worker_start"
    worker_calls = [call[0] for call in view.calls[mark:] if call[0].startswith("worker_")]
    assert worker_calls == (["worker_end_process"] if ends else []) + [new]


def test_a_build_that_fails_to_parse_abandons_the_running_process_and_a_later_one_leaves_none():
    view, worker, _, presenter = running_large_build()
    view.form = replace(big_form(55), k_text="x")

    presenter.build()

    assert (worker.ended, presenter.phase) == (1, EMPTY), "the first build was left running"
    view.form = big_form(55)
    presenter.build()
    assert (len(worker.processes), worker.ended, presenter.phase) == (2, 1, BUSY)


@pytest.mark.parametrize("close", ["quit", "shutdown"])
def test_closing_ends_the_running_process_before_the_window_closes(close):
    view, worker, _, presenter = running_large_build()

    getattr(presenter, close)()

    assert worker.ended == 1
    names = [call[0] for call in view.calls]
    if close == "quit":
        assert names.index("worker_end_process") < names.index("close")

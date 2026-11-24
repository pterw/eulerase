"""The confirmation before a large build: asked once, in the stated range; a no changes nothing."""

import pytest
from gui_fakes import SHORT_READ, FakeScheduler, FakeView, FakeWorker

from eulerase.gui.presenter import BUSY, EMPTY, READY, FormValues, Presenter

THRESHOLD = 50
K = 5


def big_form(length: int, both_strands: bool = False) -> FormValues:
    """One sequence of `length` letters: it holds `length - 4` 5-mers on one strand."""
    text = ("ACGT" * length)[:length]
    return FormValues(text, "sequences", str(K), both_strands, False, "", False)


@pytest.fixture(autouse=True)
def low_threshold(monkeypatch):
    monkeypatch.setattr("eulerase.gui.presenter.CONFIRM_ABOVE_KMERS", THRESHOLD)


def rig():
    view = FakeView(SHORT_READ)
    scheduler = FakeScheduler()
    worker = FakeWorker(view.calls)
    return view, worker, scheduler, Presenter(view, worker, scheduler)


def confirms(view):
    return [call for call in view.calls if call[0] == "confirm"]


@pytest.mark.parametrize(
    ("length", "both_strands", "kmer_count", "asks"),
    [
        (54, False, 50, False),  # equal to the threshold: no question
        (55, False, 51, True),  # one above it: asks
        (29, True, 50, False),  # 25 k-mers, counted twice, equal to the threshold
        (30, True, 52, True),  # 26 k-mers, counted twice, above it
        (1004, False, 1000, True),  # the question names the estimate, thousands separated
    ],
)
def test_a_build_asks_only_above_the_threshold_counting_both_strands_twice(
    length, both_strands, kmer_count, asks
):
    view, worker, _, presenter = rig()
    view.form = big_form(length, both_strands)
    # By hand: (length - K + 1) k-mers on a strand, doubled with both strands.
    assert kmer_count == (length - K + 1) * (2 if both_strands else 1)
    assert (kmer_count > THRESHOLD) == asks  # the cases sit on each side of the line

    presenter.build()

    assert len(confirms(view)) == int(asks)
    if asks:
        assert f"{kmer_count:,}" in confirms(view)[0][1]
    assert (len(worker.jobs), presenter.phase) == (1, BUSY)


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
    mark, jobs = len(view.calls), len(worker.jobs)

    presenter.build()

    names = [call[0] for call in view.calls[mark:]]
    assert names == ["read_form", "canvas_aspect", "confirm", "show_status"]
    assert view.status_is_error is False
    assert (presenter.phase, presenter.analysis, view.result, scheduler.pending) == before
    assert len(worker.jobs) == jobs
    assert presenter.playing == (state == "walk playing")
    if state == "build running":  # not abandoned: its outcome still arrives
        worker.finish(0, presenter)
        assert (presenter.phase, view.result is not None) == (READY, True)


def test_an_input_exactly_at_the_core_limit_is_still_confirmed(monkeypatch):
    monkeypatch.setattr("eulerase.kmers.MAX_KMERS", 51)
    view, worker, _, presenter = rig()
    view.form = big_form(55)  # 51 k-mers: above the threshold, equal to the limit

    presenter.build()

    assert len(confirms(view)) == 1
    assert (len(worker.jobs), presenter.phase) == (1, BUSY)


def test_an_input_above_the_core_limit_is_refused_not_confirmed(monkeypatch):
    monkeypatch.setattr("eulerase.kmers.MAX_KMERS", 40)
    view, worker, _, presenter = rig()
    view.form = big_form(30, True)  # 52 k-mers: above the threshold and the limit

    presenter.build()
    worker.finish(0, presenter)

    assert confirms(view) == []
    assert (presenter.phase, view.status_is_error) == (EMPTY, True)
    assert "40" in view.status

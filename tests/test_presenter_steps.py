"""Stepping through the walk and the one play timer, with a fake view, worker and scheduler."""

import dataclasses

import pytest
from gui_fakes import FakeScheduler, FakeView, FakeWorker

from eulerase.gui import presenter as presenter_module
from eulerase.gui.presenter import FormValues, Presenter

THREE_TRAILS = FormValues(
    ">a\nACGTACGGA\n>b\nTTGCCATTA\n>c\nGGGATCCC", "sequences", "4", False, False, "", False
)


class TrackingScheduler(FakeScheduler):
    """Logs every after, cancel and fire, so a test can replay what was pending when."""

    def __init__(self, events: list[tuple]):
        super().__init__()
        self.events = events

    def after(self, delay_ms, callback):
        handle = super().after(delay_ms, callback)
        self.events.append(("after", handle))
        return handle

    def cancel(self, handle):
        self.events.append(("cancel", handle))
        super().cancel(handle)

    def fire_all(self):
        self.events.extend(("fire", handle) for handle in self.pending)
        super().fire_all()


def peak_pending(events):
    """The most timers pending at once, replayed from the log alone."""
    pending, peak = set(), 0
    for kind, *rest in events:
        if kind == "after":
            pending.add(rest[0])
        elif kind in ("cancel", "fire"):
            pending.discard(rest[0])
        peak = max(peak, len(pending))
    return peak


def ready(form=THREE_TRAILS):
    """A presenter showing a finished build; returns it with its view, scheduler and layout log."""
    view = FakeView(form)
    scheduler = TrackingScheduler(view.calls)
    worker = FakeWorker(view.calls)
    layouts = []

    def layout(a, aspect):
        layouts.append(a)

    presenter = Presenter(view, worker, scheduler, layout=layout)
    presenter.build()
    worker.finish(0, presenter)
    return presenter, view, scheduler, worker, layouts


def steps_shown(view, mark=0):
    return [call[1] for call in view.calls[mark:] if call[0] == "show_step"]


def test_pressing_play_again_quickly_leaves_at_most_one_tick_pending():
    presenter, view, scheduler, *_ = ready()
    assert presenter.total > 2
    presenter.set_step(0)

    presenter.play()
    presenter.play()  # no pause between
    assert len(scheduler.pending) == 1
    presenter.pause()
    presenter.play()
    presenter.pause()
    presenter.play()
    presenter.play()
    assert (presenter.step, len(scheduler.pending)) == (0, 1)
    scheduler.fire_all()

    assert peak_pending(scheduler.events) == 1
    assert presenter.step == 1


def test_play_advances_one_step_per_tick_and_stops_at_the_end():
    presenter, view, scheduler, *_ = ready()
    presenter.set_step(0)
    mark = len(view.calls)

    presenter.play()
    for _ in range(presenter.total + 2):
        scheduler.fire_all()

    assert steps_shown(view, mark) == list(range(1, presenter.total + 1))
    assert [c for c in view.calls[mark:] if c[0] == "set_playing"] == [
        ("set_playing", True),
        ("set_playing", False),
    ]
    assert (scheduler.pending, presenter.playing) == ({}, False)


def test_play_at_the_end_restarts_from_the_first_step():
    presenter, view, scheduler, *_ = ready()
    assert presenter.step == presenter.total > 1
    mark = len(view.calls)

    presenter.play()
    assert (presenter.step, steps_shown(view, mark)) == (0, [0])
    scheduler.fire_all()

    assert presenter.step == 1


def test_set_step_clamps_to_the_walk_and_ignores_the_current_step():
    presenter, view, *_ = ready()
    total = presenter.total
    assert total > 2
    mark = len(view.calls)

    presenter.set_step(total)  # already there
    assert steps_shown(view, mark) == []
    presenter.set_step(-5)
    presenter.set_step(-5)
    presenter.set_step(999)

    assert steps_shown(view, mark) == [0, total]
    assert presenter.step == total


@pytest.mark.parametrize("action", ["build", "failed build", "shutdown"])
def test_a_new_build_a_failed_build_and_shutdown_cancel_the_timer(action):
    presenter, view, scheduler, worker, _ = ready()
    presenter.set_step(0)
    presenter.play()
    assert (len(scheduler.pending), presenter.playing) == (1, True)
    if action == "failed build":
        view.form = FormValues("ACGT", "sequences", "abc", False, False, "", False)

    {"build": presenter.build, "failed build": presenter.build, "shutdown": presenter.shutdown}[
        action
    ]()
    assert (scheduler.pending, presenter.playing) == ({}, False)  # cancelled by the action itself
    scheduler.fire_all()

    assert presenter.step == 0  # the cancelled tick never ran


def test_play_does_nothing_before_a_result_and_while_busy():
    view = FakeView(THREE_TRAILS)
    scheduler = FakeScheduler()
    worker = FakeWorker(view.calls)
    presenter = Presenter(view, worker, scheduler)

    presenter.play()  # empty
    presenter.build()
    presenter.play()  # busy: the job has started and not finished

    assert (scheduler.pending, presenter.playing) == ({}, False)
    assert ("set_playing", True) not in view.calls


def test_stepping_moves_by_the_stated_amounts_and_asks_for_no_analysis_and_no_layout():
    presenter, view, scheduler, worker, layouts = ready()
    total = presenter.total
    middle = total // 2
    assert 1 < middle < total - 1

    reached = []
    for action, start in [
        (presenter.step_back, middle),
        (presenter.step_forward, middle),
        (presenter.step_start, middle),
        (presenter.step_end, middle),
    ]:
        presenter.set_step(start)
        action()
        reached.append(presenter.step)
    assert reached == [middle - 1, middle + 1, 0, total]

    presenter.set_step(0)
    presenter.toggle_play()
    assert presenter.playing
    scheduler.fire_all()
    presenter.toggle_play()
    assert (presenter.playing, scheduler.pending) == (False, {})
    assert [c for c in view.calls if c[0] == "set_playing"][-2:] == [
        ("set_playing", True),
        ("set_playing", False),
    ]
    assert (len(worker.jobs), len(layouts)) == (1, 1)


def test_selecting_an_edge_row_jumps_to_its_global_step():
    presenter, view, *_ = ready()
    edges = next(t for t in view.result.tables if t.name == "Edges")
    local = edges.columns.index("step")
    # one global step per edge, in walk order: the per-trail step restarts, so it differs
    assert any(int(row[local]) != i + 1 for i, row in enumerate(edges.rows))

    for i in range(len(edges.rows)):
        presenter.set_step(0)
        mark = len(view.calls)
        presenter.row_selected("Edges", i)
        assert steps_shown(view, mark) == [i + 1]

    presenter.set_step(0)
    mark = len(view.calls)
    for table, index in [("Nodes", 1), ("Edges", -1), ("Edges", len(edges.rows))]:
        presenter.row_selected(table, index)
    assert (steps_shown(view, mark), presenter.step) == ([], 0)


def test_selecting_an_edge_row_reads_the_global_step_column_by_its_header(monkeypatch):
    real = presenter_module.result_view

    def reordered(a, layout):
        result = real(a, layout)
        tables = []
        for t in result.tables:
            if t.name == "Edges":  # the columns in the reverse of their usual order
                t = dataclasses.replace(
                    t, columns=t.columns[::-1], rows=tuple(r[::-1] for r in t.rows)
                )
            tables.append(t)
        return dataclasses.replace(result, tables=tuple(tables))

    monkeypatch.setattr(presenter_module, "result_view", reordered)
    presenter, view, *_ = ready()
    assert next(t for t in view.result.tables if t.name == "Edges").columns[0] == "times seen"
    presenter.set_step(0)
    mark = len(view.calls)

    presenter.row_selected("Edges", 2)

    assert steps_shown(view, mark) == [3]


def test_quit_cancels_the_timer_and_closes_the_view_once():
    presenter, view, scheduler, *_ = ready()
    presenter.set_step(0)
    presenter.play()
    handle = next(iter(scheduler.pending))

    presenter.quit()

    assert scheduler.pending == {}
    assert view.calls.count(("close",)) == 1
    assert view.calls.index(("cancel", handle)) < view.calls.index(("close",))


def test_oracle_rejects_two_pending_ticks():
    assert peak_pending([("after", 1), ("after", 2)]) == 2
    assert peak_pending([("after", 1), ("cancel", 1), ("after", 2), ("fire", 2)]) == 1

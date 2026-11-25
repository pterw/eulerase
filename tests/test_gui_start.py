"""The worker thread, the scheduler and the start function, checked without a window."""

import inspect
import threading
import time

import pytest

pytest.importorskip("tkinter", reason="the window module imports tkinter; nothing here opens one")
import tkinter  # noqa: E402

assert getattr(tkinter.Tk, "eulerase_window_guard", False) and getattr(
    tkinter.Toplevel, "eulerase_window_guard", False
), "guard not installed"

from gui_fakes import FakeScheduler, FakeView, protocol_gaps  # noqa: E402

from eulerase.errors import EuleraseError  # noqa: E402
from eulerase.gui import app  # noqa: E402
from eulerase.gui.presenter import (  # noqa: E402
    EXAMPLES,
    EXPORT_KINDS,
    Presenter,
    Scheduler,
    View,
    Worker,
)

WAIT_SECONDS = 5


def command_problems(commands) -> list[str]:
    """What is wrong with each command: no such Presenter method, or arguments it cannot take."""
    problems = []
    for c in commands:
        method = getattr(Presenter, c.method, None)
        if not callable(method):
            problems.append(f"{c.label}: no method {c.method}")
        elif not _binds(method, c.args):
            problems.append(f"{c.label}: {c.method} cannot take {c.args}")
    return problems


def _binds(method, args) -> bool:
    try:
        inspect.signature(method).bind(None, *args)
    except TypeError:
        return False
    return True


def test_every_command_names_a_method_of_the_presenter_class():
    assert command_problems(app.COMMANDS) == []
    exported = {c.args[0] for c in app.COMMANDS if c.method == "export"}
    assert exported == set(EXPORT_KINDS)
    examples = [c for c in app.COMMANDS if c.menu == "Examples"]
    assert {c.method for c in examples} == {"load_example"}
    assert [c.args for c in examples] == [(e.name,) for e in EXAMPLES]
    assert [c.method for c in app.COMMANDS if c.label == "Quit"] == ["quit"]


def test_oracle_rejects_a_command_with_no_method():
    missing = app.Command("File", "Broken", None, "no_such")
    too_many = app.Command("File", "Extra", None, "build", ("one",))
    assert len(command_problems([missing, too_many])) == 2


def wait_for_outcomes(worker, count: int = 1) -> list:
    found = []
    deadline = time.monotonic() + WAIT_SECONDS
    while len(found) < count and time.monotonic() < deadline:
        found += worker.poll()
        time.sleep(0.01)
    return found


@pytest.mark.parametrize("layout_raises", [False, True])
def test_the_worker_runs_the_job_off_the_calling_thread_and_returns_the_outcome(layout_raises):
    ran_on, daemon_flags = [], []
    marker = object()

    def layout(a, aspect):
        ran_on.append(threading.get_ident())
        daemon_flags.append(threading.current_thread().daemon)
        if layout_raises:
            raise RuntimeError("layout failed")
        return marker

    worker = app.ThreadWorker()
    Presenter(FakeView(), worker, FakeScheduler(), layout=layout).build()
    (outcome,) = wait_for_outcomes(worker)
    assert ran_on and ran_on != [threading.get_ident()]
    # a daemon thread, so closing the window during a build does not keep the process alive
    assert daemon_flags == [True]
    if layout_raises:
        assert outcome.analysis is None and "RuntimeError" in outcome.error
    else:
        assert outcome.layout is marker and outcome.error is None
    assert worker.poll() == []


class RecordingRoot:
    def __init__(self):
        self.scheduled = []

    def after(self, delay_ms, callback):
        self.scheduled.append((delay_ms, callback))
        return len(self.scheduled)


class QueuedWorker:
    def __init__(self, outcomes):
        self.outcomes = outcomes

    def poll(self):
        found, self.outcomes = self.outcomes, []
        return found


class RecordingPresenter:
    def __init__(self, fails_on=None):
        self.received = []
        self.fails_on = fails_on

    def on_build_outcome(self, outcome):
        self.received.append(outcome)
        if outcome == self.fails_on:
            raise RuntimeError("the presenter failed")


def test_the_poll_hands_each_outcome_to_the_presenter_and_polls_again():
    root, presenter = RecordingRoot(), RecordingPresenter()
    app._poll(root, QueuedWorker(["first", "second"]), presenter)
    assert presenter.received == ["first", "second"]
    ((delay, callback),) = root.scheduled
    assert delay > 0
    callback()
    assert len(root.scheduled) == 2 and presenter.received == ["first", "second"]


def test_a_failing_presenter_does_not_end_the_polling():
    root, presenter = RecordingRoot(), RecordingPresenter(fails_on="first")
    with pytest.raises(RuntimeError):
        app._poll(root, QueuedWorker(["first", "second"]), presenter)
    assert presenter.received == ["first"]
    assert len(root.scheduled) == 1, "the next look is booked before the outcomes are handed on"


def test_no_display_is_one_sentence_not_a_traceback(monkeypatch):
    def no_display(*args, **kwargs):
        raise tkinter.TclError("no display name and no $DISPLAY environment variable")

    monkeypatch.setattr(tkinter, "Tk", no_display)
    assert app.tkinter.Tk is no_display, "the fake is not in place; start is not called"
    with pytest.raises(EuleraseError) as caught:
        app.start()
    message = str(caught.value)
    assert "display" in message and "\n" not in message


def test_start_takes_one_optional_form():
    parameters = inspect.signature(app.start).parameters
    assert list(parameters) == ["form"] and parameters["form"].default is None


@pytest.mark.parametrize(
    ("cls", "protocol"),
    [(app.Window, View), (app.ThreadWorker, Worker), (app.TkScheduler, Scheduler)],
)
def test_the_window_classes_implement_the_protocols(cls, protocol):
    assert protocol_gaps(cls, protocol) == []


def test_the_scheduler_hands_timers_to_the_root():
    class Root:
        def __init__(self):
            self.cancelled = []

        def after(self, delay_ms, callback):
            return ("timer", delay_ms, callback)

        def after_cancel(self, handle):
            self.cancelled.append(handle)

    root, ticked = Root(), []
    scheduler = app.TkScheduler(root)
    handle = scheduler.after(550, ticked.append)
    assert handle == ("timer", 550, ticked.append)
    scheduler.cancel(handle)
    assert root.cancelled == [handle]

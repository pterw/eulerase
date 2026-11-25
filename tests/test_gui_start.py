"""The worker thread, the scheduler and the start function, checked without a window."""

import inspect
import multiprocessing
import threading
import time

import pytest

pytest.importorskip("tkinter", reason="the window module imports tkinter; nothing here opens one")
import tkinter  # noqa: E402

assert getattr(tkinter.Tk, "eulerase_window_guard", False) and getattr(
    tkinter.Toplevel, "eulerase_window_guard", False
), "guard not installed"

from genomes import child_env  # noqa: E402
from gui_fakes import (  # noqa: E402
    FakeScheduler,
    FakeView,
    protocol_gaps,
    sized_job,
    sleeping_job,
    small_job,
    unreadable_job,
    vanishing_job,
)

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
PROCESS_WAIT_SECONDS = 30  # a spawned child starts a fresh interpreter and imports gui_fakes
SLEEP_SECONDS = 3


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


def wait_for_outcomes(worker, count: int = 1, seconds: float = WAIT_SECONDS) -> list:
    found = []
    deadline = time.monotonic() + seconds
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


@pytest.fixture
def process_worker(tmp_path, monkeypatch):
    """A worker whose spawned child takes the suite's child environment and works in tmp_path;
    whatever it left running is ended, and every child is gone when the test is over."""
    for name, value in child_env().items():
        monkeypatch.setenv(name, value)
    monkeypatch.chdir(tmp_path)
    worker = app.ThreadWorker()
    yield worker
    worker.end_process()
    for child in multiprocessing.active_children():  # a finished build's child is reaped late
        child.join(PROCESS_WAIT_SECONDS)


@pytest.mark.parametrize("child_gone_first", [False, True])
def test_a_process_build_returns_its_outcome_with_its_token_once(process_worker, child_gone_first):
    text = "x" * 1_000_000  # far more than a start-up pipe holds
    # the interpreter's first start pays a one-time cost (the resource tracker), so time the second
    process_worker.start_process(6, small_job, ("first",))
    assert [o.token for o in wait_for_outcomes(process_worker, seconds=PROCESS_WAIT_SECONDS)] == [6]
    for first_child in multiprocessing.active_children():
        first_child.join(PROCESS_WAIT_SECONDS)
    started = time.perf_counter()
    process_worker.start_process(7, sized_job, (text,))
    call_seconds = time.perf_counter() - started
    (child,) = multiprocessing.active_children()
    if child_gone_first:  # sent and exited before the first look: still a result, not a loss
        child.join(PROCESS_WAIT_SECONDS)
        assert not child.is_alive()

    (outcome,) = wait_for_outcomes(process_worker, seconds=PROCESS_WAIT_SECONDS)
    total_seconds = time.perf_counter() - started

    assert (outcome.token, outcome.error) == (7, str(len(text)))
    assert process_worker.poll() == []
    # started with `spawn` (never `fork` beside Tk), and as a daemon
    assert isinstance(child, multiprocessing.get_context("spawn").Process) and child.daemon
    # the window's thread does not wait for the child to start and read its argument
    assert total_seconds > 0.2, "no real child start-up to compare with"
    assert call_seconds < total_seconds / 4


def test_starting_a_process_ends_the_one_still_running(process_worker, tmp_path):
    process_worker.start_process(1, sleeping_job, (str(tmp_path / "started"), SLEEP_SECONDS))
    (first,) = multiprocessing.active_children()

    process_worker.start_process(2, small_job, ("second",))

    assert not first.is_alive(), "the first build was left running"
    (outcome,) = wait_for_outcomes(process_worker, seconds=PROCESS_WAIT_SECONDS)
    assert (outcome.token, outcome.error) == (2, "second")


def test_ending_a_process_leaves_no_child_and_never_delivers_its_outcome(process_worker, tmp_path):
    marker = tmp_path / "started"
    process_worker.start_process(3, sleeping_job, (str(marker), SLEEP_SECONDS))
    deadline = time.monotonic() + PROCESS_WAIT_SECONDS
    while not marker.exists() and time.monotonic() < deadline:
        time.sleep(0.01)
    seen = time.monotonic()
    assert marker.exists() and multiprocessing.active_children(), "the job has not started"

    process_worker.end_process()

    deadline = seen + SLEEP_SECONDS / 2
    while multiprocessing.active_children() and time.monotonic() < deadline:
        time.sleep(0.01)
    # gone well before the job's own sleep would have ended it: it was ended, not left to finish
    assert multiprocessing.active_children() == []
    time.sleep(max(0, seen + SLEEP_SECONDS + 0.5 - time.monotonic()))  # past that sleep's end
    assert process_worker.poll() == []
    process_worker.end_process()  # a second end, and one with nothing started, do nothing
    app.ThreadWorker().end_process()


@pytest.mark.parametrize("job", [vanishing_job, unreadable_job])  # no outcome; one unreadable
def test_a_process_with_no_readable_outcome_gives_one_failure_with_its_token(process_worker, job):
    process_worker.start_process(5, job, ())

    (outcome,) = wait_for_outcomes(process_worker, seconds=PROCESS_WAIT_SECONDS)

    assert (outcome.token, outcome.analysis, outcome.layout) == (5, None, None)
    assert outcome.error and "\n" not in outcome.error and "Traceback" not in outcome.error
    assert process_worker.poll() == []


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

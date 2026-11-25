"""Fakes for the presenter's three collaborators, and a check that they match the protocols.

Imported flat, as `genomes.py` is. Nothing here touches a toolkit.
"""

import inspect
import os
import time
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path

from eulerase.analysis import Analysis
from eulerase.gui.presenter import BuildOutcome, FormValues, ResultView

SHORT_READ = FormValues("ACGTACGTTG", "sequences", "3", False, False, "", False)


class FakeView:
    """Records every call, in order, as a tuple of the method name and its arguments."""

    def __init__(self, form: FormValues = SHORT_READ, aspect: float = 2.0):
        self.calls: list[tuple] = []
        self.form = form
        self.aspect = aspect
        self.result: ResultView | None = None
        self.busy = False
        self.status: str | None = None
        self.status_is_error = False
        self.exports: Mapping[str, bool] = {}
        self.k_enabled: bool | None = None
        self.confirm_answer = True
        self.open_path: str | None = None
        self.save_path: str | None = None
        self.export_answer = True

    def _record(self, *call) -> None:
        self.calls.append(call)

    def read_form(self) -> FormValues:
        self._record("read_form")
        return self.form

    def canvas_aspect(self) -> float:
        self._record("canvas_aspect")
        return self.aspect

    def set_form(self, form: FormValues) -> None:
        self._record("set_form", form)
        self.form = form

    def set_k_enabled(self, enabled: bool) -> None:
        self._record("set_k_enabled", enabled)
        self.k_enabled = enabled

    def show_status(self, message: str, error: bool = False) -> None:
        self._record("show_status", message, error)
        self.status, self.status_is_error = message, error

    def show_busy(self, busy: bool) -> None:
        self._record("show_busy", busy)
        self.busy = busy

    def show_result(self, result: ResultView) -> None:
        self._record("show_result", result)
        self.result = result

    def clear_result(self) -> None:
        self._record("clear_result")
        self.result = None

    def show_step(self, step: int, total: int) -> None:
        self._record("show_step", step, total)

    def set_playing(self, playing: bool) -> None:
        self._record("set_playing", playing)

    def set_exports_enabled(self, enabled: Mapping[str, bool]) -> None:
        self._record("set_exports_enabled", dict(enabled))
        self.exports = dict(enabled)

    def confirm(self, message: str) -> bool:
        self._record("confirm", message)
        return self.confirm_answer

    def ask_open_path(self) -> str | None:
        self._record("ask_open_path")
        return self.open_path

    def ask_save_path(self, kind: str, suggested: str) -> str | None:
        self._record("ask_save_path", kind, suggested)
        return self.save_path

    def confirm_export(self, paths: Sequence[str], existing: Sequence[str]) -> bool:
        self._record("confirm_export", tuple(paths), tuple(existing))
        return self.export_answer

    def save_figure(self, analysis: Analysis, path: str) -> None:
        self._record("save_figure", analysis, path)

    def close(self) -> None:
        self._record("close")


class FakeScheduler:
    def __init__(self):
        self.pending: dict[int, Callable[[], None]] = {}
        self._next = 0

    def after(self, delay_ms: int, callback: Callable[[], None]) -> int:
        self._next += 1
        self.pending[self._next] = callback
        return self._next

    def cancel(self, handle: object) -> None:
        self.pending.pop(handle, None)

    def fire_all(self) -> None:
        due, self.pending = list(self.pending.values()), {}
        for callback in due:
            callback()


class FakeWorker:
    """Keeps the started jobs and processes; nothing runs until a test says so. A shared `log`
    records the starts and the ends in order with a view's calls."""

    def __init__(self, log: list[tuple] | None = None):
        self.jobs: list[Callable[[], BuildOutcome]] = []
        self.processes: list[tuple[int, Callable[..., BuildOutcome], tuple]] = []
        self.ended = 0  # the calls of end_process
        self._log = log

    def start(self, job: Callable[[], BuildOutcome]) -> None:
        if self._log is not None:
            self._log.append(("worker_start",))
        self.jobs.append(job)

    def start_process(self, token: int, function: Callable[..., BuildOutcome], args: tuple) -> None:
        if self._log is not None:
            self._log.append(("worker_start_process",))
        self.processes.append((token, function, args))

    def end_process(self) -> None:
        if self._log is not None:
            self._log.append(("worker_end_process",))
        self.ended += 1

    def run_without_delivering(self, i: int) -> BuildOutcome:
        return self.jobs[i]()

    def finish(self, i: int, presenter) -> None:
        presenter.on_build_outcome(self.run_without_delivering(i))


def small_job(token: int, text: str) -> BuildOutcome:
    """A finished build for the real-process tests; `text` rides in `error` so it can be compared.
    Module-level, because a spawned child imports it by name."""
    return BuildOutcome(token, None, None, text)


def sized_job(token: int, text: str) -> BuildOutcome:
    """Takes an argument of any size and answers with its length, so the outcome stays small."""
    return BuildOutcome(token, None, None, str(len(text)))


def sleeping_job(token: int, marker: str, seconds: float) -> BuildOutcome:
    """Writes `marker` when it starts, then sleeps: a build still running."""
    Path(marker).write_text("started")
    time.sleep(seconds)
    return BuildOutcome(token, None, None, "slept")


def vanishing_job(token: int) -> BuildOutcome:
    """Ends its process at once, with no outcome."""
    os._exit(3)


def fail_to_rebuild():
    raise ValueError("this object cannot be rebuilt")


class Unreadable:
    """Pickles in the child; unpickling it in the parent raises an error that is not a lost
    pipe."""

    def __reduce__(self):
        return (fail_to_rebuild, ())


def unreadable_job(token: int) -> Unreadable:
    """Finishes, but what it sends cannot be read back."""
    return Unreadable()


def protocol_gaps(cls: type, protocol: type) -> list[str]:
    """Each public method of the protocol that `cls` lacks, or whose parameter names differ."""
    gaps = []
    for name, wanted in vars(protocol).items():
        if name.startswith("_") or not inspect.isfunction(wanted):
            continue
        found = getattr(cls, name, None)
        if not callable(found):
            gaps.append(f"{name}: missing")
        elif list(inspect.signature(found).parameters) != list(
            inspect.signature(wanted).parameters
        ):
            gaps.append(f"{name}: parameter names differ")
    return gaps

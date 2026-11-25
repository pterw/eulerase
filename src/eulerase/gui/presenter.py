"""What the window decides, without a toolkit: shapes, protocols and the pure functions."""

import os
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Protocol

from eulerase import analysis, io, kmers, sequence
from eulerase.analysis import Analysis, Settings
from eulerase.errors import EuleraseError
from eulerase.sequence import Record

EMPTY = "empty"
BUSY = "busy"
READY = "ready"
EXPORT_KINDS = ("figure", "csv", "fasta", "genbank", "gfa")
CONFIRM_ABOVE_KMERS = 100_000
PLAY_INTERVAL_MS = 550
DEFAULT_WINDOW = (1280, 640)
MIN_WINDOW = (900, 560)
# Tk's geometry counts the client area only: the title bar, menu bar, borders and task bar
# take this much of the screen on top of it.
SCREEN_ALLOWANCE = (0, 120)

_TEXT_MEASURES = ("settings", "orf rule")
_CONTIG_EXPORTS = ("fasta", "genbank", "gfa")
# kind: (the dialog's suggested name, the writer of its text; None where the window draws it)
_TEXT_EXPORTS = {
    "figure": ("eulerase.png", None),
    "fasta": ("contigs.fasta", io.fasta_text),
    "genbank": ("contigs.gbk", io.genbank_text),
    "gfa": ("assembly.gfa", io.gfa_text),
}


@dataclass(frozen=True)
class FormValues:
    text: str
    input_mode: str  # "sequences" or "kmers"
    k_text: str  # exactly what the user typed
    both_strands: bool
    unique: bool
    min_multiplicity_text: str  # exactly what the user typed; empty means none
    clip_tips: bool


@dataclass(frozen=True)
class TableView:
    name: str
    columns: tuple[str, ...]
    rows: tuple[tuple[str, ...], ...]
    note: str | None


@dataclass(frozen=True)
class ResultView:
    analysis: Analysis
    layout: object
    headline: str
    status: str
    details: tuple[str, ...]
    tables: tuple[TableView, ...]
    steps: int


@dataclass(frozen=True)
class BuildOutcome:
    """Exactly one of `analysis` and `error` is not None."""

    token: int
    analysis: Analysis | None
    layout: object
    error: str | None


@dataclass(frozen=True)
class Example:
    name: str
    text: str
    input_mode: str
    k: int
    both_strands: bool
    unique: bool
    min_multiplicity: int | None
    clip_tips: bool


def _example(name: str, text: str, k: int, **settings) -> Example:
    fields = {
        "input_mode": "sequences",
        "both_strands": False,
        "unique": False,
        "min_multiplicity": None,
        "clip_tips": False,
    }
    return Example(name, text, k=k, **(fields | settings))


EXAMPLES = (
    _example("Lemon text, k=3 (single path)", "LEMONLEMONADE", 3),
    _example("DNA with a repeat, k=4", "ACGTTGCAACGTTGCA", 4),
    _example("Tandem repeat (self-loop), k=3", "AAAAAC", 3),
    _example("Two unrelated reads (2 trails), k=4", ">read1\nACGTACGGA\n>read2\nTTGCCATTA", 4),
    _example("Both strands, k=5", "ATGGCGTACGTTAGC", 5, both_strands=True),
    _example("A k-mer list that closes into a circuit", "ACG CGT GTA TAC", 3, input_mode="kmers"),
    _example(
        "Overlapping reads (assembly), k=5",
        ">r1\nACGTGCAATG\n>r2\nCAATGCCGTA\n>r3\nCCGTAACGTG",
        5,
    ),
    _example(
        "Reads with one error, tip clipping on, k=4",
        ">ok\nACGTGCGTACGT\n>err\nACGTTT",
        4,
        clip_tips=True,
    ),
)


class View(Protocol):
    def read_form(self) -> FormValues: ...

    def canvas_aspect(self) -> float:
        """Canvas width over height, read on the window's thread; the figure aspect when
        the width or the height is under 2 pixels (a canvas not yet mapped)."""
        ...

    def set_form(self, form: FormValues) -> None:
        """Put the values in the controls; build nothing."""
        ...

    def set_k_enabled(self, enabled: bool) -> None: ...

    def show_status(self, message: str, error: bool = False) -> None: ...

    def show_busy(self, busy: bool) -> None: ...

    def show_result(self, result: ResultView) -> None: ...

    def clear_result(self) -> None: ...

    def show_step(self, step: int, total: int) -> None:
        """The first `step` steps of the walk are traversed."""
        ...

    def set_playing(self, playing: bool) -> None: ...

    def set_exports_enabled(self, enabled: Mapping[str, bool]) -> None:
        """Keys are EXPORT_KINDS."""
        ...

    def confirm(self, message: str) -> bool: ...

    def ask_open_path(self) -> str | None: ...

    def ask_save_path(self, kind: str, suggested: str) -> str | None: ...

    def confirm_export(self, paths: Sequence[str], existing: Sequence[str]) -> bool: ...

    def save_figure(self, analysis: Analysis, path: str) -> None:
        """May raise EuleraseError (a wrong extension) or OSError (an unwritable path or a
        missing directory: render.save_figure does not turn those into EuleraseError), so
        the caller catches both."""
        ...

    def close(self) -> None:
        """End the window: the one way the presenter closes it."""
        ...


class Worker(Protocol):
    def start(self, job: Callable[[], BuildOutcome]) -> None:
        """Run the job off the window's thread; its outcome later goes to
        Presenter.on_build_outcome on the window's thread."""
        ...


class Scheduler(Protocol):
    def after(self, delay_ms: int, callback: Callable[[], None]) -> object: ...

    def cancel(self, handle: object) -> None: ...


def _table_rows(a: Analysis, name: str) -> list[tuple[str, ...]]:
    rows = analysis.rows(a, name)
    if name == "Assembly":  # the numeric rows only; the text rows go to the details
        return [row for row in rows if row[0] not in _TEXT_MEASURES]
    return rows


def _orf_rule_line(a: Analysis) -> str:
    return f"orf rule: {dict(analysis.rows(a, 'Assembly'))['orf rule']}"


def _table_view(a: Analysis, name: str) -> TableView:
    shown, note = analysis.limit_rows(_table_rows(a, name))
    if name == "Contigs":
        note = "\n".join(part for part in (note, _orf_rule_line(a)) if part)
    return TableView(name, analysis.columns(name), tuple(shown), note)


def table_views(a: Analysis) -> tuple[TableView, ...]:
    return tuple(_table_view(a, name) for name in analysis.TABLES)


def result_view(a: Analysis, layout: object) -> ResultView:
    settings = analysis.settings_text(a.settings)
    return ResultView(
        analysis=a,
        layout=layout,
        headline=a.headline,
        status=f"{a.headline} | {settings}",
        details=(f"settings: {settings}", _orf_rule_line(a), *a.notes),
        tables=table_views(a),
        steps=sum(len(trail) for trail in a.trails),
    )


def records_to_text(records: Sequence[Record]) -> str:
    return "\n".join(f">{r.name}\n{r.sequence}" for r in records)


def estimate_kmers(
    records: Sequence[Record], k: int | None, *, both_strands: bool, input_mode: str
) -> int:
    """The k-mers a build would generate, counted before any duplicate is dropped."""
    if input_mode == "kmers":
        count = len(records)
    else:
        count = sum(max(0, len(r.sequence) - k + 1) for r in records)
    return count * (2 if both_strands else 1)


def fit_window(screen_width: int, screen_height: int) -> tuple[int, int]:
    """The default window cut to the screen less the allowance, never below the minimum."""
    return tuple(
        max(least, min(default, screen - taken))
        for least, default, screen, taken in zip(
            MIN_WINDOW, DEFAULT_WINDOW, (screen_width, screen_height), SCREEN_ALLOWANCE, strict=True
        )
    )


def _whole_number(text: str, what: str) -> int:
    try:
        return int(text)  # int() already ignores spaces around the digits
    except ValueError:
        raise EuleraseError(
            f"{what} must be a whole number, not {text.strip()!r}; type digits only."
        ) from None


def _settings(form: FormValues) -> Settings:
    minimum = form.min_multiplicity_text.strip()
    common = {
        "both_strands": form.both_strands,
        "unique": form.unique,
        "min_multiplicity": _whole_number(minimum, "The minimum multiplicity") if minimum else None,
        "clip_tips": form.clip_tips,
        "input_mode": form.input_mode,
    }
    if form.input_mode == "kmers":  # k comes from the k-mers; the box is never read
        return Settings(k=None, **common)
    return Settings(k=kmers.validate_k(_whole_number(form.k_text, "k")), **common)


def _unwritable(error: OSError, after: str = "") -> str:
    reason = error.strerror or "unwritable"
    named = f" {error.filename}" if error.filename else ""
    return (
        f"The file{named} could not be written ({reason}){after}; check the folder and permissions."
    )


def _records(form: FormValues) -> list[Record]:
    parse = sequence.parse_kmer_list if form.input_mode == "kmers" else sequence.parse_records
    return parse(form.text)


class Presenter:
    """The window's decisions. It reaches the window only through `View`."""

    def __init__(
        self,
        view: View,
        worker: Worker,
        scheduler: Scheduler,
        *,
        layout: Callable[[Analysis, float], object] | None = None,
    ):
        self._view = view
        self._worker = worker
        self._scheduler = scheduler
        self._layout = layout or (lambda a, aspect: None)
        self._phase = EMPTY
        self._analysis: Analysis | None = None
        self._step = 0
        self._total = 0
        self._token = 0
        self._tick: object | None = None  # the one pending play timer
        self._result: ResultView | None = None  # what the window shows

    @property
    def phase(self) -> str:
        return self._phase

    @property
    def analysis(self) -> Analysis | None:
        return self._analysis

    @property
    def step(self) -> int:
        return self._step

    @property
    def total(self) -> int:
        return self._total

    @property
    def playing(self) -> bool:
        return self._tick is not None

    def mode_changed(self, input_mode: str) -> None:
        self._view.set_k_enabled(input_mode == "sequences")

    def load_example(self, name: str) -> None:
        example = next(e for e in EXAMPLES if e.name == name)
        minimum = example.min_multiplicity
        self.open_with(
            FormValues(
                example.text,
                example.input_mode,
                str(example.k),
                example.both_strands,
                example.unique,
                "" if minimum is None else str(minimum),
                example.clip_tips,
            )
        )

    def open_with(self, form: FormValues | None) -> None:
        if form is None:
            self.load_example(EXAMPLES[0].name)
            return
        self._view.set_form(form)
        self.mode_changed(form.input_mode)
        self.build()

    def build(self) -> None:
        form = self._view.read_form()
        aspect = self._view.canvas_aspect()
        try:
            settings, records = _settings(form), _records(form)
        except EuleraseError as error:
            self._fail(str(error))
            return
        if not self._confirmed(records, settings):
            self._view.show_status("Build cancelled; nothing was changed.")
            return
        self.pause()
        self._token += 1  # abandons any build still running
        self._clear_result()
        self._phase = BUSY
        self._view.show_busy(True)
        self._worker.start(self._job(records, settings, aspect, self._token))

    def _confirmed(self, records: Sequence[Record], settings: Settings) -> bool:
        """Ask once about a large build; above the core's limit the core refuses, so no question."""
        count = estimate_kmers(
            records,
            settings.k,
            both_strands=settings.both_strands,
            input_mode=settings.input_mode,
        )
        if not CONFIRM_ABOVE_KMERS < count <= kmers.MAX_KMERS:
            return True
        return self._view.confirm(
            f"This build makes {count:,} k-mers and may take a while. Build anyway?"
        )

    def stop(self) -> None:
        # No tick to cancel: BUSY is entered only by build(), which pauses, and play() does
        # nothing outside READY.
        if self._phase != BUSY:
            return
        self._token += 1
        self._phase = EMPTY
        self._view.show_busy(False)
        self._view.show_status("Stopped.")

    def on_build_outcome(self, outcome: BuildOutcome) -> None:
        if outcome.token != self._token or self._phase != BUSY:
            return
        if outcome.error is not None:
            self._fail(outcome.error)
            return
        result = result_view(outcome.analysis, outcome.layout)
        self._analysis = outcome.analysis
        self._phase = READY
        self._step = self._total = result.steps
        self._result = result
        self._view.show_result(result)
        self._view.show_step(self._total, self._total)
        self._view.show_busy(False)
        self._view.set_playing(False)
        self._view.set_exports_enabled({kind: self._enabled(kind) for kind in EXPORT_KINDS})
        self._view.show_status(result.status)

    def open_file(self) -> None:
        path = self._view.ask_open_path()
        if path:
            self.load_file(path)

    def load_file(self, path: str) -> None:
        try:
            records = io.read_records(path)
        except EuleraseError as error:  # the screen stays as it was
            self._view.show_status(str(error), error=True)
            return
        loaded = replace(
            self._view.read_form(), text=records_to_text(records), input_mode="sequences"
        )
        self.open_with(loaded)

    def _enabled(self, kind: str) -> bool:
        """Whether `kind` can be exported from the analysis shown: a build has succeeded."""
        return kind not in _CONTIG_EXPORTS or bool(self._analysis.contigs)

    def export(self, kind: str) -> None:
        if self._phase != READY or not self._enabled(kind):
            return
        try:
            message = self._export_csv() if kind == "csv" else self._export_one(kind)
        except EuleraseError as error:
            self._view.show_status(str(error), error=True)
        except OSError as error:
            self._view.show_status(_unwritable(error), error=True)
        else:
            if message:
                self._view.show_status(message)

    def _export_one(self, kind: str) -> str | None:
        suggested, writer = _TEXT_EXPORTS[kind]
        path = self._view.ask_save_path(kind, suggested)
        if not path:
            return None
        if writer is None:
            self._view.save_figure(self._analysis, path)
        else:
            io.write_text_file(path, writer(self._analysis))
        return f"Saved {path}."

    def _export_csv(self) -> str | None:
        chosen = self._view.ask_save_path("csv", "eulerase.csv")
        if not chosen:
            return None
        prefix = chosen[:-4] if chosen.lower().endswith(".csv") else chosen
        paths = list(io.csv_paths(prefix).values())
        existing = [path for path in paths if os.path.exists(path)]
        if not self._view.confirm_export(paths, existing):
            return None
        written = 0
        try:
            for table, path in zip(analysis.TABLES, paths, strict=True):
                io.write_text_file(path, io.csv_text(self._analysis, table))
                written += 1
        except OSError as error:
            after = f", after {written} of {len(paths)} files were written"
            raise EuleraseError(_unwritable(error, after)) from error
        return f"Saved {len(paths)} files: {prefix}_*.csv."

    def _job(
        self, records: list[Record], settings: Settings, aspect: float, token: int
    ) -> Callable[[], BuildOutcome]:
        layout = self._layout

        def job() -> BuildOutcome:
            try:
                a = analysis.analyze(records, settings)
                return BuildOutcome(token, a, layout(a, aspect), None)
            except EuleraseError as error:
                return BuildOutcome(token, None, None, str(error))
            except Exception as error:  # the busy state must always end
                message = (
                    f"The build failed with a {type(error).__name__}: this is a bug in Eulerase."
                )
                return BuildOutcome(token, None, None, message)

        return job

    def set_step(self, n: int) -> None:
        n = max(0, min(n, self._total))
        if n != self._step:
            self._step = n
            self._view.show_step(n, self._total)

    def step_back(self) -> None:
        self.set_step(self._step - 1)

    def step_forward(self) -> None:
        self.set_step(self._step + 1)

    def step_start(self) -> None:
        self.set_step(0)

    def step_end(self) -> None:
        self.set_step(self._total)

    def play(self) -> None:
        if self._phase != READY:
            return
        self._cancel_tick()  # at most one tick is ever pending
        if self._step == self._total:
            self.set_step(0)
        self._view.set_playing(True)
        self._schedule_tick()

    def pause(self) -> None:
        was_playing = self.playing
        self._cancel_tick()
        if was_playing:
            self._view.set_playing(False)

    def toggle_play(self) -> None:
        if self.playing:
            self.pause()
        else:
            self.play()

    def row_selected(self, table: str, index: int) -> None:
        if self._result is None or table != "Edges":
            return
        edges = next(t for t in self._result.tables if t.name == "Edges")
        if 0 <= index < len(edges.rows):
            self.set_step(int(edges.rows[index][edges.columns.index("global step")]))

    def shutdown(self) -> None:
        self._cancel_tick()

    def quit(self) -> None:
        self.shutdown()
        self._view.close()

    def _schedule_tick(self) -> None:
        self._tick = self._scheduler.after(PLAY_INTERVAL_MS, self._on_tick)

    def _on_tick(self) -> None:
        self._tick = None
        self.set_step(self._step + 1)
        if self._step < self._total:
            self._schedule_tick()
        else:
            self._view.set_playing(False)

    def _cancel_tick(self) -> None:
        if self._tick is not None:
            self._scheduler.cancel(self._tick)
            self._tick = None

    def _clear_result(self) -> None:
        self._analysis = None
        self._step = self._total = 0
        self._result = None
        self._view.clear_result()
        self._view.set_exports_enabled(dict.fromkeys(EXPORT_KINDS, False))

    def _fail(self, message: str) -> None:
        self.pause()
        self._phase = EMPTY
        self._clear_result()
        self._view.show_busy(False)
        self._view.show_status(message, error=True)

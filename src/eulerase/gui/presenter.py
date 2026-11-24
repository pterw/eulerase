"""What the window decides, without a toolkit: shapes, protocols and the pure functions."""

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Protocol

from eulerase import analysis
from eulerase.analysis import Analysis
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

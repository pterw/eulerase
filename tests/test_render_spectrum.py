"""Tests for the k-mer spectrum chart: one bar per multiplicity, the minimum marked on it."""

from collections import Counter

import pytest
from genomes import brute_kmers, random_dna
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_rgb
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from eulerase import render
from eulerase.analysis import Settings, analyze
from eulerase.sequence import Record

K = 3
READS = tuple(random_dna(70, seed) for seed in (1, 2))
BOTH = pytest.mark.parametrize("both", [False, True], ids=["one-strand", "both-strands"])


def expected_spectrum(both: bool) -> Counter:
    """Multiplicity to number of distinct k-mers, counted from the reads alone."""
    return Counter(Counter(brute_kmers(READS, K, both)).values())


def chart(both: bool, minimum: int | None = None, reads=READS, figsize=(6.4, 4.8)):
    records = [Record(f"r{i}", s) for i, s in enumerate(reads)]
    a = analyze(records, Settings(k=K, both_strands=both, min_multiplicity=minimum))
    fig = Figure(figsize=figsize)
    FigureCanvasAgg(fig)
    render.draw_spectrum(fig, a)
    return a, fig, fig.axes[0]


def shown_ticks(ax) -> list[tuple[object, Text]]:
    """(axis, label) of every tick that lies in the view, after a draw."""
    ax.figure.canvas.draw()
    return [
        (axis, label)
        for axis in (ax.xaxis, ax.yaxis)
        for loc, label in zip(axis.get_majorticklocs(), axis.get_majorticklabels())
        if min(axis.get_view_interval()) <= loc <= max(axis.get_view_interval())
    ]


def bars(ax) -> list[tuple[int, int, Rectangle]]:
    """Each bar as (multiplicity, height, patch), read from the axes."""
    return [
        (round(p.get_x() + p.get_width() / 2), round(p.get_height()), p)
        for p in ax.patches
        if isinstance(p, Rectangle)
    ]


def drawn_spectrum(ax) -> dict[int, int]:
    return {x: h for x, h, _ in bars(ax)}


@BOTH
def test_bars_show_the_spectrum_of_the_input(both):
    want = expected_spectrum(both)
    # two multiplicities with different heights, so a swap of x and height would show
    assert len(set(want.values())) > 1 and len(want) > 1
    _, _, ax = chart(both)
    assert drawn_spectrum(ax) == dict(want)


def test_oracle_rejects_a_wrong_bar():
    want = expected_spectrum(False)
    _, _, ax = chart(False)
    ax.patches[0].set_height(ax.patches[0].get_height() + 1)
    assert drawn_spectrum(ax) != dict(want)


@BOTH
def test_the_minimum_is_marked_and_the_bars_it_removes_differ(both):
    minimum = 3
    counts = Counter(brute_kmers(READS, K, both))
    removed = sum(1 for c in counts.values() if c < minimum)
    assert 0 < removed < len(counts)  # the chart has both kinds of bar
    a, _, ax = chart(both, minimum)
    assert a.removed_kmers == removed
    cut = [(x, h, p) for x, h, p in bars(ax) if x < minimum]
    kept = [(x, h, p) for x, h, p in bars(ax) if x >= minimum]
    assert cut and kept
    assert sum(h for _, h, _ in cut) == removed
    cut_style = {(to_rgb(p.get_facecolor()), p.get_hatch()) for _, _, p in cut}
    kept_style = {(to_rgb(p.get_facecolor()), p.get_hatch()) for _, _, p in kept}
    assert len(cut_style) == 1 and len(kept_style) == 1
    (cut_color, cut_hatch), (kept_color, kept_hatch) = cut_style.pop(), kept_style.pop()
    assert cut_color != kept_color and cut_hatch and not kept_hatch
    assert (cut_color, kept_color) == (to_rgb(render.FUTURE), to_rgb(render.TRAIL_COLORS[0]))
    lines = [line.get_xdata() for line in ax.lines if line.get_xdata()[0] == line.get_xdata()[-1]]
    assert any(minimum - 1 < xs[0] <= minimum for xs in lines)  # at the minimum, or half a bar left
    assert [str(minimum) in t.get_text() for t in ax.texts] == [True]  # the line is labelled
    fig = ax.figure
    render.draw_spectrum(fig, a)  # a redraw on the same figure replaces the chart
    assert len(fig.axes) == 1


@pytest.mark.parametrize(
    "reads, both, bars_drawn",
    [(READS, False, len(expected_spectrum(False))), (("ACGTT",), False, 1)],
    ids=["many-bars", "one-bar"],
)
def test_ticks_are_whole_numbers(reads, both, bars_drawn):
    _, _, ax = chart(both, reads=reads)
    assert len(drawn_spectrum(ax)) == bars_drawn  # one bar is the case that gave fractional ticks
    ticks = shown_ticks(ax)
    assert {axis for axis, _ in ticks} == {ax.xaxis, ax.yaxis}  # both axes were read
    assert all(label.get_text().isdigit() for _, label in ticks)


@pytest.mark.parametrize("minimum", [None, 50])
def test_every_text_lies_inside_a_small_figure(minimum):
    if minimum is not None:
        assert minimum > max(expected_spectrum(False))  # above every multiplicity
    _, fig, ax = chart(False, minimum, figsize=(4.5, 3))
    shown = [ax.xaxis.label, ax.yaxis.label, *ax.texts, *(t for _, t in shown_ticks(ax))]
    assert len(ax.texts) == (minimum is not None)
    for t in shown:
        box = t.get_window_extent()
        assert fig.bbox.x0 <= box.x0 and box.x1 <= fig.bbox.x1, t.get_text()
        assert fig.bbox.y0 <= box.y0 and box.y1 <= fig.bbox.y1, t.get_text()


def test_no_mark_one_colour_and_no_hatch_without_a_minimum():
    _, _, ax = chart(False)
    assert not ax.lines
    assert len({to_rgb(p.get_facecolor()) for _, _, p in bars(ax)}) == 1
    assert not any(p.get_hatch() for _, _, p in bars(ax))


@pytest.mark.parametrize("minimum", [None, 3])
def test_chart_text_names_no_coverage_or_depth(minimum):
    _, fig, ax = chart(False, minimum)
    texts = [t.get_text() for t in fig.findobj(lambda o: hasattr(o, "get_text"))]
    texts += [ax.get_xlabel(), ax.get_ylabel(), ax.get_title()]
    joined = " ".join(texts).lower()
    assert "times seen" in joined  # the text was read: the words below are really absent
    assert not any(word in joined for word in ("coverage", "depth", "cov="))

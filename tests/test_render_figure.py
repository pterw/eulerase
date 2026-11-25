"""Tests for the redrawn figure: copy marks, the current step, trail colours, the legend, the frame.

The new names are reached through `render.` so that a missing name fails the test that uses it.
"""

import re
from collections import Counter

import pytest
from figure_probe import (
    INPUTS,
    SAVED_SIZE,
    SMALL_SIZE,
    built,
    dash,
    drawn,
    elements,
    extent,
    frame,
    inside,
    legends,
    overlap,
    ten,
    visible_texts,
)
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_rgb
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from eulerase import render
from eulerase.analysis import Settings, analyze, settings_text
from eulerase.sequence import Record

NARROW_SIZE = (4.5, 3.6)  # inches: narrower than the headline and the settings line at one line
SIZES = {"saved": SAVED_SIZE, "small": SMALL_SIZE, "narrow": NARROW_SIZE}


# Inputs with a long settings line (a minimum multiplicity is set), for the frame and the legend:
# the reads and the k of each.
LONG_SETTINGS = {
    "one": (("LEMONLEMONADE",), 3),
    "three": (("ACGTACGGA", "TTGCCATTAGGC", "GGATCCTTAAGC"), 5),
}


@pytest.fixture
def input_named(monkeypatch):
    """The analysis of an input by name; the ten-trail input has more nodes than the overview
    limit, which is raised here so that it is drawn in full."""

    def make(name):
        if name == "ten":
            monkeypatch.setattr(render, "OVERVIEW_ABOVE_NODES", 1000)
            return ten()
        if name in LONG_SETTINGS:
            return _long_settings(*LONG_SETTINGS[name])
        return INPUTS[name]()

    return make


def _steps(a):
    return [s for trail in a.trails for s in trail]


def _copies(a, upto=None):
    """How many times each (trail, k-mer) is met in the first `upto` steps (all by default)."""
    return Counter((s.trail, s.kmer) for s in _steps(a) if upto is None or s.global_step <= upto)


def _by_gid(fig, role):
    return {artist.get_gid(): artist for artist in elements(fig, role)}


def _rgb(artist):
    return artist.get_facecolor()[:3]


def _edge(artist):
    return artist.get_edgecolor()[:3]


def _long_settings(reads, k):
    """An analysis whose settings line is long, with nothing filtered out by it (a minimum of 1)."""
    records = [Record(f"r{i}", text) for i, text in enumerate(reads, 1)]
    a = analyze(records, Settings(k=k, min_multiplicity=1))
    assert len(settings_text(a.settings)) > 60
    return a


# ---------------------------------------------------------------------------
# Copy marks, the count, the current step, the colour of each arrow
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ("repeated", "shared"))
def test_a_repeated_kmer_is_one_arrow_with_one_mark_per_copy_and_the_used_copies_are_filled(name):
    a = INPUTS[name]()
    total = _copies(a)
    if name == "repeated":
        assert max(total.values()) == 3
    else:  # a k-mer seen more than once in two trails: keyed by k-mer alone, they would merge
        repeated = Counter(kmer for (_, kmer), n in total.items() if n > 1)
        assert max(repeated.values()) >= 2
    for upto in range(len(_steps(a)) + 1):
        fig = built(a, upto=upto)
        used = _copies(a, upto)
        assert sorted(_by_gid(fig, "arrow")) == sorted(f"arrow:{t}:{k}" for t, k in total)
        marks = Counter()
        filled = Counter()
        for mark in elements(fig, "mark"):
            marks[mark.get_gid()] += 1
            filled[mark.get_gid()] += mark.get_fill()
        for (trail, kmer), copies in total.items():
            gid = f"mark:{trail}:{kmer}"
            # a k-mer seen once has no line of marks: its arrow already says whether it is used
            assert marks[gid] == (copies if copies > 1 else 0), (upto, gid)
            assert filled[gid] == (used[(trail, kmer)] if copies > 1 else 0), (upto, gid)


@pytest.mark.parametrize("copies", (render.MAX_COPY_MARKS, render.MAX_COPY_MARKS + 1))
def test_above_the_copy_limit_the_marks_become_a_count(copies):
    text = "".join("ABC" + letter for letter in "XYZUVWPQR"[:copies])
    a = analyze([Record("r", text)], Settings(k=3))
    assert _copies(a)[(1, "ABC")] == copies
    for upto in range(len(_steps(a)) + 1):
        fig = built(a, upto=upto)
        marks = [m for m in elements(fig, "mark") if m.get_gid() == "mark:1:ABC"]
        counts = elements(fig, "count")
        if copies <= render.MAX_COPY_MARKS:
            assert len(marks) == copies and counts == []
        else:
            used = _copies(a, upto)[(1, "ABC")]
            assert marks == []
            (count,) = counts
            assert count.get_gid() == "count:1:ABC"
            assert [int(n) for n in re.findall(r"\d+", count.get_text())] == [used, copies]


@pytest.mark.parametrize("copies", (render.MAX_COPY_MARKS, 100))
def test_the_marks_and_the_count_stay_inside_their_labels_slot(copies, monkeypatch):
    # a two-letter k-mer has a narrow label: its marks at full pitch, or its count at full size,
    # would be wider than the slot the label reserved, so they must shrink to fit
    reserved = []
    place = render.place_labels
    monkeypatch.setattr(
        render, "place_labels", lambda *args: reserved.append(place(*args)) or reserved[-1]
    )
    a = analyze([Record("r", "AB" * copies)], Settings(k=2))
    fig, _ = drawn(a, SMALL_SIZE)
    (placement,) = reserved
    axes = fig.axes[0]
    to_px = fig.dpi / 72
    scale = axes.get_window_extent().height / to_px  # the plot's height in points
    slots = {}
    for key, placed in placement.placed.items():
        low, high = (axes.transData.transform((x / scale, 0))[0] for x in placed.box[::2])
        slots[key.removeprefix("label:")] = (low, high)
    inked = elements(fig, "mark") + elements(fig, "count")
    assert {artist.get_gid().split(":")[0] for artist in inked} == {
        "mark" if copies <= render.MAX_COPY_MARKS else "count"
    }
    first, last = slots["1:AB"]
    if copies <= render.MAX_COPY_MARKS:
        unshrunk = copies * render.MARK_PITCH * to_px
    else:
        unshrunk = len(f"{copies}/{copies}") * render.MONO_ADVANCE * render.COUNT_PT * to_px
    assert unshrunk > last - first
    for artist in inked:
        low, high = slots[artist.get_gid().split(":", 1)[1]]
        box = extent(fig, artist)
        assert low - 0.5 <= box[0] and box[2] <= high + 0.5, artist.get_gid()


@pytest.mark.parametrize("name", ("lemon", "repeated", "shared"))
def test_the_current_kmer_and_node_are_marked(name):
    a = INPUTS[name]()
    steps = _steps(a)
    for upto in range(len(steps) + 1):
        fig = built(a, upto=upto)
        bands = _by_gid(fig, "band")
        assert len(bands) == len(elements(fig, "arrow"))  # all exist from the start
        shown = {gid: b for gid, b in bands.items() if b.get_visible()}
        nodes = _by_gid(fig, "node")
        lit = [gid for gid, n in nodes.items() if _rgb(n) == to_rgb(render.LEMON)]
        labels = _by_gid(fig, "label")
        lemon_labels = [
            gid for gid, t in labels.items() if _rgb(t.get_bbox_patch()) == to_rgb(render.LEMON)
        ]
        if upto == 0:
            assert (list(shown), lit, lemon_labels) == ([], [], [])
            continue
        step = steps[upto - 1]
        assert step.global_step == upto
        assert list(shown) == [f"band:{step.trail}:{step.kmer}"]
        assert _edge(shown[f"band:{step.trail}:{step.kmer}"]) == to_rgb(render.LEMON)
        assert lit == [f"node:{step.target}"]
        assert lemon_labels == [f"label:{step.trail}:{step.kmer}"]
        outline = nodes[f"node:{step.target}"].get_linewidth()
        assert outline > max(n.get_linewidth() for gid, n in nodes.items() if gid not in lit)


@pytest.mark.parametrize("name", ("lemon", "shared", "ten"))
def test_a_used_arrow_takes_its_trails_colour_and_dash_and_an_unused_one_is_grey(name, input_named):
    a = input_named(name)
    assert len(a.trails) == {"lemon": 1, "shared": 3, "ten": 10}[name]
    total = _copies(a)
    steps = range(len(_steps(a)) + 1)
    if name == "ten":  # every trail used, one step at a time would be 135 figures
        steps = sorted({0, len(_steps(a)), *(trail[0].global_step for trail in a.trails)})
    for upto in steps:
        arrows = _by_gid(built(a, upto=upto), "arrow")
        used = _copies(a, upto)
        thin, thick = [], []
        for trail, kmer in total:
            arrow = arrows[f"arrow:{trail}:{kmer}"]
            if used[(trail, kmer)]:
                assert _edge(arrow) == to_rgb(render.TRAIL_COLORS[(trail - 1) % 3]), (upto, trail)
                assert dash(arrow) == render.TRAIL_DASHES[((trail - 1) // 3) % 3], (upto, trail)
                thick.append(arrow.get_linewidth())
            else:
                assert _edge(arrow) == to_rgb(render.FUTURE), (upto, trail, kmer)
                thin.append(arrow.get_linewidth())
        if thin and thick:
            assert min(thick) > max(thin)


# ---------------------------------------------------------------------------
# The legend
# ---------------------------------------------------------------------------

# The colour and dash of trails 1 to 10, typed out: the dash changes at trails 4 and 7, and trail
# 10 has the colour and the dash of trail 1.
TEN_TRAILS = [
    ("#2A78D6", "solid"),
    ("#EB6834", "solid"),
    ("#1BAF7A", "solid"),
    ("#2A78D6", "dashed"),
    ("#EB6834", "dashed"),
    ("#1BAF7A", "dashed"),
    ("#2A78D6", "dotted"),
    ("#EB6834", "dotted"),
    ("#1BAF7A", "dotted"),
    ("#2A78D6", "solid"),
]


def test_legend_entries_follow_the_trails_and_the_colour_rule():
    a = ten()
    entries = render.legend_entries(a)
    trails = [e for e in entries if e.kind == "trail"]
    assert [e.kind for e in entries] == ["trail"] * 10 + ["start", "end"]
    for number, (entry, (colour, line)) in enumerate(zip(trails, TEN_TRAILS), 1):
        assert re.search(rf"\b{number}\b", entry.label), entry.label
        assert (entry.color, entry.dash) == (colour, line), number
    start, end = entries[10:]
    assert (start.color, start.dash) == ("#1F1F1E", "solid")
    assert (end.color, end.dash) == ("#1F1F1E", "dashed")


@pytest.mark.parametrize("name", ("lemon", "shared"))
def test_the_figure_legend_shows_the_entries(name):
    a = INPUTS[name]()
    entries = render.legend_entries(a)
    (legend,) = legends(built(a))
    assert [t.get_text() for t in legend.get_texts()] == [e.label for e in entries]
    handles = legend.legend_handles
    assert len(handles) == len(entries) > 2
    for handle, entry in zip(handles, entries):
        assert to_rgb(handle.get_color()) == to_rgb(entry.color)
        assert dash(handle) == entry.dash


def legend_defects(fig):
    """What is wrong with the figure's legend: none, several, one over the plot or over a text."""
    found = legends(fig)
    if len(found) != 1:
        return ["no legend" if not found else "several legends"]
    box = extent(fig, found[0])
    defects = ["over the plot"] if overlap(box, extent(fig, fig.axes[0]), margin=0) else []
    return defects + [
        ("over text", t.get_text())
        for t in visible_texts(fig)
        if overlap(box, extent(fig, t), margin=0)
    ]


@pytest.mark.parametrize("size", SIZES.values(), ids=SIZES.keys())
@pytest.mark.parametrize("name", ("one", "three"))
def test_the_legend_lies_outside_the_plot_and_touches_no_text(name, size, input_named):
    fig, _ = drawn(input_named(name), size)
    assert legend_defects(fig) == []


def _forged(legend_at=None, headline_at=0.9):
    """A figure with a plot axes, a headline and, when asked for, a legend at `legend_at`."""
    fig = Figure(figsize=(4, 3))
    FigureCanvasAgg(fig)
    fig.add_axes((0.1, 0.3, 0.8, 0.5))
    fig.text(0.5, headline_at, "A headline", ha="center", va="top", gid="headline")
    if legend_at is not None:
        sample = [Line2D([0], [0], label="trail 1")]
        fig.legend(handles=sample, loc="center", bbox_to_anchor=legend_at)
    fig.canvas.draw()
    return fig


def test_oracle_rejects_an_overlapping_legend():
    assert legend_defects(_forged((0.5, 0.05))) == []
    assert legend_defects(_forged()) == ["no legend"]
    assert legend_defects(_forged((0.5, 0.55))) == ["over the plot"]
    over_headline = legend_defects(_forged((0.5, 0.9)))
    assert over_headline == [("over text", "A headline")]


# ---------------------------------------------------------------------------
# The frame: headline, reading line, settings line, left-out text, plot area
# ---------------------------------------------------------------------------


def frame_defects(fig):
    """What is wrong with the figure's frame: a text or the legend outside the figure, a missing
    text, or two texts that touch."""
    boxes = {}
    found = []
    for role in ("headline", "reading", "settings"):
        texts = frame(fig, role)
        if len(texts) != 1:
            found.append((role, "not one text"))
        boxes.update({f"{role}{i}": extent(fig, t) for i, t in enumerate(texts)})
    boxes.update({f"left-out{i}": extent(fig, t) for i, t in enumerate(frame(fig, "left-out"))})
    boxes.update({"legend": extent(fig, legend) for legend in legends(fig)})
    found += [(name, "outside") for name, box in boxes.items() if not inside(box, fig)]
    plot = extent(fig, fig.axes[0])
    found += [
        (name, "plot")
        for name, box in boxes.items()
        if name.startswith("left-out") and overlap(box, plot, margin=0)
    ]
    names = sorted(boxes)
    for i, first in enumerate(names):
        found += [
            (first, other)
            for other in names[i + 1 :]
            if overlap(boxes[first], boxes[other], margin=0)
        ]
    return found


@pytest.mark.parametrize("size", SIZES.values(), ids=SIZES.keys())
@pytest.mark.parametrize("name", ("one", "three", "ten"))
def test_the_frame_lies_inside_the_figure(name, size, input_named):
    fig, _ = drawn(input_named(name), size)
    assert legends(fig)
    assert frame_defects(fig) == []


def test_oracle_rejects_a_frame_outside_the_figure():
    missing = [("reading", "not one text"), ("settings", "not one text")]
    assert frame_defects(_forged((0.5, 0.05))) == missing
    fig = _forged((0.5, 0.05))
    for role, at in (("reading", 0.8), ("settings", 0.75)):
        fig.text(0.5, at, role, ha="center", va="top", gid=role)
    fig.canvas.draw()
    assert frame_defects(fig) == []
    fig.texts[0].set_text("A headline that is far wider than the figure that holds it " * 2)
    fig.canvas.draw()
    assert ("headline0", "outside") in frame_defects(fig)
    fig.texts[0].set_text("A headline")
    fig.texts[1].set_position((0.5, 0.9))  # on the headline
    fig.canvas.draw()
    assert ("headline0", "reading0") in frame_defects(fig)
    fig.texts[1].set_position((0.5, 0.8))
    fig.text(0.5, 0.55, "Labels: 1 left out.", ha="center", gid="left-out")  # over the plot
    fig.canvas.draw()
    assert frame_defects(fig) == [("left-out0", "plot")]


def _shortened_only():
    """Three copies of one k-mer, and a canvas (the small size) where its marks line fits nowhere
    but its text does, and no label is left out."""
    return analyze([Record("r", "CAABACAACA")], Settings(k=3))


LEFT_OUT_CASES = {  # name: (analysis, size, (some label shortened, some label left out))
    "none": (INPUTS["lemon"], SAVED_SIZE, (False, False)),
    "both": (INPUTS["overlap"], SMALL_SIZE, (True, True)),
    "only-left-out": (INPUTS["three"], SMALL_SIZE, (False, True)),
    "both-shared": (INPUTS["shared"], SMALL_SIZE, (True, True)),
    "only-shortened": (_shortened_only, SMALL_SIZE, (True, False)),
}


@pytest.mark.parametrize("case", LEFT_OUT_CASES)
def test_the_left_out_text_appears_only_when_a_count_is_above_zero(case):
    make, size, shape = LEFT_OUT_CASES[case]
    a = make()
    fig, _ = drawn(a, size)
    total = _copies(a)
    wanted = len(total) + 2 * len(a.trails)  # a label for each arrow, a start and an end word
    labels = _by_gid(fig, "label")
    left_out = wanted - len(labels) - len(elements(fig, "word"))
    drawn_marks = Counter(m.get_gid().removeprefix("mark:") for m in elements(fig, "mark"))
    counted = {g.removeprefix("count:") for g in _by_gid(fig, "count")}
    shortened = sum(
        1
        for (trail, kmer), copies in total.items()
        if copies > 1
        and f"label:{trail}:{kmer}" in labels
        and f"{trail}:{kmer}" not in drawn_marks
        and f"{trail}:{kmer}" not in counted
    )
    # the shape of the case, counted from the figure: a text shown for the wrong reason or a count
    # swapped with the other must show
    assert (shortened > 0, left_out > 0) == shape
    assert shortened != left_out or not shortened
    expected = [[shortened, left_out]] if any(shape) else []
    shown = [[int(n) for n in re.findall(r"\d+", t.get_text())] for t in frame(fig, "left-out")]
    assert shown == expected


def _plot_aspect(fig):
    box = extent(fig, fig.axes[0])
    return (box[2] - box[0]) / (box[3] - box[1])


def aspect_defects(fig):
    """What is wrong with the plot area's aspect: it is not FIGURE_ASPECT within one percent."""
    return [] if _plot_aspect(fig) == pytest.approx(render.FIGURE_ASPECT, rel=0.01) else ["aspect"]


@pytest.mark.parametrize("size", SIZES.values(), ids=SIZES.keys())
def test_the_plot_area_of_a_saved_figure_has_the_figure_aspect(size):
    fig, _ = drawn(INPUTS["shared"](), size)
    assert aspect_defects(fig) == []


def test_oracle_rejects_a_plot_area_of_the_wrong_aspect():
    assert aspect_defects(_forged()) == ["aspect"]  # 0.8 by 0.5 of a 4 by 3 inch figure
    square = Figure(figsize=(4, 4))
    FigureCanvasAgg(square)
    square.add_axes((0.1, 0.1, 0.8, 0.8))
    square.canvas.draw()
    assert aspect_defects(square) == ["aspect"]

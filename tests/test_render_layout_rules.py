"""Rules of the figure's plot: where labels, rings and words go, and the overview."""

import dataclasses

import pytest
from figure_probe import INPUTS, SAVED_SIZE, SMALL_SIZE, drawn, elements, extent, inside, overlap
from genomes import random_dna
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_rgb
from matplotlib.figure import Figure
from matplotlib.patches import Circle

from eulerase import render
from eulerase.analysis import Settings, analyze
from eulerase.sequence import Record

SIZES = {"saved": SAVED_SIZE, "small": SMALL_SIZE}


def _texts(fig):
    """The visible labels and words of the plot, by gid."""
    return {
        artist.get_gid(): artist
        for role in ("label", "word")
        for artist in elements(fig, role)
        if artist.get_visible()
    }


def collisions(fig):
    """The defects of the texts of a figure: one on a node, on another text, or outside."""
    boxes = {gid: extent(fig, artist) for gid, artist in _texts(fig).items()}
    discs = {disc.get_gid(): extent(fig, disc) for disc in elements(fig, "node")}
    found = []
    for gid, box in boxes.items():
        found += [(gid, node) for node, disc in discs.items() if overlap(box, disc)]
        if not inside(box, fig):
            found.append((gid, "outside"))
    gids = sorted(boxes)
    for i, first in enumerate(gids):
        found += [(first, g) for g in gids[i + 1 :] if overlap(boxes[first], boxes[g])]
    return found


def positions(fig):
    """The position and visibility of every label and word, by gid."""
    return {
        artist.get_gid(): (artist.get_position(), artist.get_visible())
        for role in ("label", "word")
        for artist in elements(fig, role)
    }


def moves(first, second):
    """The gids whose position or visibility differ between two figures."""
    before, after = positions(first), positions(second)
    return sorted(gid for gid in before.keys() | after.keys() if before.get(gid) != after.get(gid))


def ring_defects(fig, a, layout):
    """The defects of each trail's rings and of the words that were placed (a word that fits
    nowhere is left out and counted, the ring stays) against the layout points of its ends."""
    found = []
    for t, trail in enumerate(a.trails, 1):
        for kind, node, style in (
            ("start", trail[0].source, ("-", "solid")),
            ("end", trail[-1].target, ("--", "dashed")),
        ):
            (ring,) = [e for e in elements(fig, "ring") if e.get_gid() == f"ring:{kind}:{t}"]
            words = [e for e in elements(fig, "word") if e.get_gid() == f"word:{kind}:{t}"]
            if ring.center != pytest.approx(layout[node]):
                found.append((kind, t, "centre"))
            if ring.get_linestyle() not in style:
                found.append((kind, t, "line style"))
            if words and words[0].get_text().split() != [kind, str(t)]:
                found.append((kind, t, "word"))
    rings = [(ring.center, ring.radius) for ring in elements(fig, "ring")]
    found += [
        ("rings", centre, "same radius")
        for centre, radius in rings
        if rings.count((centre, radius)) > 1
    ]
    return found


@pytest.mark.parametrize("size", SIZES.values(), ids=SIZES.keys())
@pytest.mark.parametrize("name", INPUTS)
def test_no_label_touches_a_node_or_another_label(name, size):
    fig, _ = drawn(INPUTS[name](), size)
    assert elements(fig, "label")
    assert collisions(fig) == []


def test_no_label_touches_a_node_or_another_label_at_any_step():
    a = INPUTS["lemon"]()
    for upto in range(sum(len(trail) for trail in a.trails) + 1):
        fig, _ = drawn(a, SAVED_SIZE, upto)
        assert elements(fig, "label"), upto
        assert collisions(fig) == [], upto


def test_labels_do_not_move_between_steps():
    a = INPUTS["lemon"]()
    first, _ = drawn(a, SAVED_SIZE, 0)
    assert len(positions(first)) > 1
    for upto in range(1, sum(len(trail) for trail in a.trails) + 1):
        later, _ = drawn(a, SAVED_SIZE, upto)
        assert moves(first, later) == [], upto


@pytest.mark.parametrize("size", SIZES.values(), ids=SIZES.keys())
@pytest.mark.parametrize("name", INPUTS)
def test_every_trail_has_a_start_ring_an_end_ring_and_a_word(name, size):
    a = INPUTS[name]()
    assert all((t[0].source == t[-1].target) == (name == "circuit") for t in a.trails)
    if name == "three":
        assert len({len(trail) for trail in a.trails}) == len(a.trails) == 3
    fig, layout = drawn(a, size)
    assert ring_defects(fig, a, layout) == []
    assert len(elements(fig, "ring")) == 2 * len(a.trails)
    if size == SAVED_SIZE:  # room for every word; on the small canvas some are left out
        assert len(elements(fig, "word")) == 2 * len(a.trails)


def _pair_count(a):
    return len({(s.source, s.target) for trail in a.trails for s in trail})


def test_above_the_overview_limit_the_plot_has_no_text_and_one_arrow_per_pair_of_nodes(
    monkeypatch,
):
    dna = random_dna(200, seed=1)  # a second read repeats every pair of nodes in another trail
    a = analyze([Record("s", dna), Record("t", dna)], Settings(k=5))
    assert len(a.trails) > 1
    steps = sum(len(trail) for trail in a.trails)
    nodes = len(a.nodes)
    assert nodes > render.OVERVIEW_ABOVE_NODES
    assert _pair_count(a) != steps  # a repeated k-mer: counting steps instead of pairs fails

    monkeypatch.setattr(render, "OVERVIEW_ABOVE_NODES", nodes)
    full, _ = drawn(a)
    assert full.axes[0].texts

    monkeypatch.setattr(render, "OVERVIEW_ABOVE_NODES", nodes - 1)
    overview, _ = drawn(a)
    assert list(overview.axes[0].texts) == []
    assert len(elements(overview, "arrow")) == _pair_count(a)
    assert len(elements(overview, "node")) == nodes
    at_start, _ = drawn(a, upto=0)  # no walk yet: every dot is ink, not the surface of a full node
    assert {to_rgb(n.get_facecolor()) for n in elements(at_start, "node")} == {to_rgb(render.INK)}


def _plot_of(a, size=SMALL_SIZE):
    """The figure and the plot that `_build_plot` returns, for the plot area's aspect."""
    fig = Figure(figsize=size)
    plot = render._build_plot(
        fig.add_subplot(111), a, render.compute_layout(a), render.FIGURE_ASPECT
    )
    return fig, plot


def test_the_left_out_count_is_the_texts_missing_from_the_plot():
    a = analyze([Record("s", random_dna(30, seed=3))], Settings(k=8))
    fig, plot = _plot_of(a)
    wanted = len(elements(fig, "arrow")) + 2 * len(
        a.trails
    )  # a label each, a start and an end word
    shown = len(elements(fig, "label")) + len(elements(fig, "word"))
    assert plot.placement.left_out == wanted - shown > 0
    assert plot.placement.left_out != plot.placement.shortened  # the two counts cannot be swapped


def _copies(a):
    """Copies of each (trail, k-mer): an arrow stands for all of them."""
    seen = {}
    for trail in a.trails:
        for s in trail:
            seen[f"label:{s.trail}:{s.kmer}"] = seen.get(f"label:{s.trail}:{s.kmer}", 0) + 1
    return seen


def test_a_label_is_shortened_to_its_k_mer_only_when_its_marks_line_does_not_fit():
    a = INPUTS["repeat"]()
    copies = _copies(a)
    assert max(copies.values()) > 1 and min(copies.values()) == 1
    _, roomy = _plot_of(a, SAVED_SIZE)
    _, tight = _plot_of(a, (4.0, 2.4))
    for plot in (roomy, tight):
        placed = plot.placement.placed
        labels = {gid: p for gid, p in placed.items() if gid.startswith("label:")}
        assert not any(p.marks for gid, p in labels.items() if copies[gid] == 1)
        cut = [gid for gid, p in labels.items() if copies[gid] > 1 and not p.marks]
        assert plot.placement.shortened == len(cut)
    assert roomy.placement.shortened == 0 < tight.placement.shortened


def _forged(label_at):
    """A figure whose plot holds a disc at the middle and one label at `label_at` (data units)."""
    fig = Figure(figsize=SAVED_SIZE)
    FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    disc = Circle((50, 50), 5)
    disc.set_gid("node:AB")
    ax.add_patch(disc)
    ax.text(*label_at, "ACG").set_gid("label:1:ACG")
    fig.canvas.draw()
    return fig


def test_oracle_rejects_a_label_on_a_node():
    assert collisions(_forged((20, 20))) == []
    assert collisions(_forged((50, 50))) == [("label:1:ACG", "node:AB")]


def test_oracle_rejects_a_label_that_moves():
    assert moves(_forged((20, 20)), _forged((20, 20))) == []
    assert moves(_forged((20, 20)), _forged((30, 20))) == ["label:1:ACG"]


def test_oracle_rejects_a_ring_off_its_node():
    a = INPUTS["lemon"]()
    fig, layout = drawn(a)
    assert ring_defects(fig, a, layout) == []
    start = a.trails[0][0].source
    off = dict(layout, **{start: (layout[start][0] + 0.2, layout[start][1])})
    assert ring_defects(fig, a, off) == [("start", 1, "centre")]


# place_labels with sizes that make every box easy to compute by hand: a character is 10 wide
# and a line 10 high, and the padding is 3, so the 8-character text fills a box 86 wide and
# 16 high around its centre, and 16 + MARKS_LINE_H high when its marks line is kept. A spot is
# (x, y, normal x, normal y, gap): the box goes beside the anchor, on the side of the normal.
BOUNDS = (0.0, 0.0, 400.0, 200.0)
KMER = "ABCDEFGH"
FIRST = (100.0, 50.0, 0.0, 1.0, 0.0)  # the box goes above (100, 50): x 57 to 143, y 50 to 66+
SECOND = (300.0, 50.0, 0.0, 1.0, 0.0)
THIN_HIT = (60.0, 67.0, 140.0, 80.0)  # meets the box with a marks line, not the one without
FULL_HIT = (60.0, 55.0, 140.0, 80.0)  # meets both
TOP = (100.0, 183.0, 0.0, 1.0, 0.0)  # the box without marks ends at 199, with them past 200


@pytest.mark.parametrize(
    "blocked, labels, kept, shortened, left_out",
    [
        ([], [([FIRST], True)], {0: True}, 0, 0),
        # a k-mer seen once has no marks line to lose: never shortened, whatever is near
        ([THIN_HIT], [([FIRST], False)], {0: False}, 0, 0),
        ([THIN_HIT], [([FIRST], True)], {0: False}, 1, 0),
        ([FULL_HIT], [([FIRST], True)], {0: None}, 0, 1),
        ([FULL_HIT], [([FIRST, SECOND], True)], {0: True}, 0, 0),
        ([], [([TOP], True)], {0: False}, 1, 0),
        ([], [([(10.0, 50.0, 0.0, 1.0, 0.0)], True)], {0: None}, 0, 1),
        # the two counts differ and each counts its own labels
        (
            [THIN_HIT],
            [([FIRST], True), ([(10.0, 50.0, 0.0, 1.0, 0.0)], True)],
            {0: False, 1: None},
            1,
            1,
        ),
    ],
    ids=[
        "free",
        "seen once",
        "shortened",
        "left out",
        "first free spot",
        "bounds shorten",
        "bounds leave out",
        "both counts",
    ],
)
def test_a_label_that_does_not_fit_loses_its_marks_line_then_is_left_out_and_counted(
    blocked, labels, kept, shortened, left_out
):
    assert render.MARKS_LINE_H > 1  # the boxes with and without a marks line differ
    wanted = [
        render.Wanted(
            f"label:{i}:{KMER}", KMER, tuple(render.Spot(*s) for s in spots), 10.0, 10.0, marks
        )
        for i, (spots, marks) in enumerate(labels)
    ]
    placement = render.place_labels(wanted, blocked, BOUNDS)
    assert (placement.shortened, placement.left_out) == (shortened, left_out)
    for i, marks in kept.items():
        got = placement.placed.get(wanted[i].key)
        assert (got and got.marks) == marks
        assert (got is None) == (marks is None)


def test_a_word_never_pushes_a_label_out_of_its_place(monkeypatch):
    """The labels of the k-mers are placed before the words of the rings: with every word left
    out (its ring still takes room), each label is where it was."""
    a = INPUTS["overlap"]()
    _, with_words = _plot_of(a, SAVED_SIZE)
    assert any(gid.startswith("word:") for gid in with_words.placement.placed)
    real = render._word
    monkeypatch.setattr(render, "_word", lambda *args: dataclasses.replace(real(*args), spots=()))
    _, without_words = _plot_of(a, SAVED_SIZE)
    assert not any(gid.startswith("word:") for gid in without_words.placement.placed)

    def label_boxes(plot):
        return {g: p.box for g, p in plot.placement.placed.items() if g.startswith("label:")}

    assert label_boxes(with_words) and label_boxes(with_words) == label_boxes(without_words)

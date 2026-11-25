"""Tests for `StepScene`: the plot drawn once and restyled for each step of the walk."""

import math
import re

import pytest
from figure_probe import (
    INPUTS,
    SAVED_SIZE,
    SMALL_SIZE,
    collisions,
    drawn,
    elements,
    extent,
    frame,
    inside,
    positions,
)
from genomes import random_dna
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.patches import FancyArrowPatch, Patch
from matplotlib.text import Text

from eulerase import render
from eulerase.analysis import Settings, analyze
from eulerase.sequence import Record

CANVASES = {"default": (988, 314), "minimum": (652, 244)}  # pixels of the window's graph canvas


def _one(text, k):
    return analyze([Record("r", text)], Settings(k=k))


def _many():
    a = _one("ABC" * 7, 3)  # the 3-mer ABC is seen seven times: above the copy marks
    assert max(sum(1 for s in t if s.kmer == "ABC") for t in a.trails) > render.MAX_COPY_MARKS
    return a


def _overview():
    a = _one(random_dna(150, seed=7), 6)
    assert len(a.nodes) > render.OVERVIEW_ABOVE_NODES
    return a


SCENE_INPUTS = {
    "lemon": INPUTS["lemon"],
    "three": INPUTS["three"],
    "many": _many,
    "overview": _overview,
}


def _total(a):
    return sum(len(trail) for trail in a.trails)


def _points(total):
    """Every step of a short walk; about ten, with both ends, of a long one."""
    return sorted({0, total, *range(0, total + 1, max(1, total // 10))})


def _scene(a, layout, size, scene_class=render.StepScene):
    """A scene on an Agg figure of `size` inches."""
    fig = Figure(figsize=size)
    FigureCanvasAgg(fig)
    return scene_class(fig, a, layout), fig


def _scene_like(fresh, a, layout, scene_class=render.StepScene):
    """A scene on a figure of the size of the plot area of `fresh`, a figure of `draw_figure`."""
    box = fresh.axes[0].get_position()
    size = (box.width * fresh.get_figwidth(), box.height * fresh.get_figheight())
    return _scene(a, layout, size, scene_class)


def _r(value):
    return round(float(value), 3)


def _style(art):
    """The properties that differ between a used, a current and a not yet used element."""
    style = {"zorder": art.get_zorder(), "alpha": art.get_alpha()}
    if isinstance(art, Text):
        style |= {"text": art.get_text(), "colour": to_rgba(art.get_color())}
        style |= {"weight": art.get_fontweight()}
        if (box := art.get_bbox_patch()) is not None:
            style["box"] = tuple(map(_r, box.get_facecolor()))
    if isinstance(art, Patch):
        style |= {"edge": tuple(map(_r, art.get_edgecolor()))}
        style |= {"face": tuple(map(_r, art.get_facecolor())), "fill": art.get_fill()}
        style |= {"dash": art.get_linestyle(), "width": _r(art.get_linewidth())}
    if isinstance(art, FancyArrowPatch):
        style["head"] = _r(art.get_mutation_scale())
    return tuple(sorted((key, repr(value)) for key, value in style.items()))


def _where(art):
    if isinstance(art, FancyArrowPatch):
        (x0, y0), (x1, y1) = art._posA_posB  # the ends, in data units, before any drawing
        return tuple(_r(v) for v in (x0, y0, x1, y1, art.get_connectionstyle().rad))
    if isinstance(art, Text):
        return tuple(_r(v) for v in art.get_position())
    return tuple(_r(v) for v in (*art.center, art.radius))


def snapshot(fig):
    """Every visible line, patch and text of the plot axes with a gid: what it is, where it is and
    its style. A hidden element is left out: the scene makes from the start, hidden, what a
    drawing makes only from some step on."""
    axes = fig.axes[0]
    found = [
        (art.get_gid(), type(art).__name__, _where(art), _style(art))
        for art in axes.get_children()
        if isinstance(art, (Patch, Text)) and art.get_gid() and art.get_visible()
    ]
    return sorted(found, key=repr)


def _fresh_and_scene(a, n, scene_class=render.StepScene, size=SMALL_SIZE):
    """The snapshots of `draw_figure` at step `n` and of a scene shown at step `n`."""
    layout = render.compute_layout(a, render.FIGURE_ASPECT)
    fresh, _ = drawn(a, size, upto=n, layout=layout)
    scene, fig = _scene_like(fresh, a, layout, scene_class)
    scene.show(n)
    assert fig.get_figwidth() / fig.get_figheight() == pytest.approx(render.FIGURE_ASPECT)
    return snapshot(fresh), snapshot(fig)


@pytest.mark.parametrize("name", SCENE_INPUTS)
def test_a_step_creates_no_figure_element(name):
    a = SCENE_INPUTS[name]()
    scene, fig = _scene(a, render.compute_layout(a, 1.5), (6, 4))
    before = fig.findobj()  # held, so that no new element can take the identity of a freed one
    nodes = [row.node for row in a.nodes]
    steps = list(range(_total(a) + 1))
    for n in steps + steps[::-1]:
        scene.show(n)
    for node in (*nodes[:3], None):
        scene.show_hover(node)
    assert {id(art) for art in fig.findobj()} == {id(art) for art in before}


@pytest.mark.parametrize("name", ("lemon", "overview"))
def test_a_step_never_calls_the_layout(name, monkeypatch):
    a = SCENE_INPUTS[name]()
    layout = render.compute_layout(a, 1.5)

    def refuse(*args, **kwargs):
        raise AssertionError("the scene asked for a layout")

    monkeypatch.setattr(render, "compute_layout", refuse)
    scene, _ = _scene(a, layout, (6, 4))
    for n in range(_total(a) + 1):
        scene.show(n)
    scene.node_at(100, 100)
    scene.show_hover(a.nodes[0].node)
    scene.show_hover(None)


@pytest.mark.parametrize("name", ("lemon", "three"))
def test_a_step_moves_no_label(name):
    a = SCENE_INPUTS[name]()
    scene, fig = _scene(a, render.compute_layout(a, 1.5), (6, 4))
    first = positions(fig)
    assert len(first) > 1 and all(visible for _, visible in first.values())
    for n in range(1, _total(a) + 1):
        scene.show(n)
        assert positions(fig) == first, n


@pytest.mark.parametrize("name", SCENE_INPUTS)
def test_a_scene_at_each_step_matches_a_fresh_drawing(name):
    a = SCENE_INPUTS[name]()
    layout = render.compute_layout(a, render.FIGURE_ASPECT)
    fresh, _ = drawn(a, SMALL_SIZE, layout=layout)
    scene, fig = _scene_like(fresh, a, layout)
    points = _points(_total(a))
    fresh_at = {n: snapshot(drawn(a, SMALL_SIZE, upto=n, layout=layout)[0]) for n in points}
    assert points[0] == 0
    assert snapshot(fig) == fresh_at[0]  # styled at once: the window shows it before any step
    for n in points + points[::-1]:
        scene.show(n)
        assert snapshot(fig) == fresh_at[n], n


def test_oracle_rejects_a_scene_that_ignores_the_step():
    a = SCENE_INPUTS["lemon"]()
    last = _total(a)

    class Ignoring(render.StepScene):
        def show(self, upto):
            pass

    first, _ = _fresh_and_scene(a, 0)
    final, _ = _fresh_and_scene(a, last)
    assert first != final  # a drawing at the first step is not one at the last
    assert _fresh_and_scene(a, last, Ignoring)[1] != final


def test_a_scene_with_no_layout_says_why(monkeypatch):
    a = SCENE_INPUTS["lemon"]()
    monkeypatch.setattr(render, "MAX_DRAW_NODES", 2)
    assert render.compute_layout(a) is None
    fig = Figure(figsize=SMALL_SIZE)
    render.draw_figure(fig, a, None)
    scene, plain = _scene(a, None, SMALL_SIZE)
    wanted = [t.get_text() for t in elements(fig, "note")]
    assert len(wanted) == 1 and "over the limit of 2" in wanted[0]
    (shown,) = [t.get_text() for t in elements(plain, "note")]
    # the window has tables and points at them; the saved figure has none and does not
    assert "Use the tables below" in shown and "Use the tables below" not in wanted[0]
    # the window has no header, so its note carries the headline; the saved figure's header does
    assert a.headline in shown and a.headline not in wanted[0]
    assert wanted[0].removesuffix(")") in shown
    # no line of the note leaves the small canvas
    assert inside(extent(plain, elements(plain, "note")[0]), plain)
    scene.show(3)
    scene.show_hover(a.nodes[0].node)
    scene.show_hover(None)
    assert scene.left_out == (0, 0) and scene.node_at(10, 10) is None


@pytest.mark.parametrize("name", INPUTS)
@pytest.mark.parametrize("canvas", CANVASES.values(), ids=CANVASES.keys())
def test_the_scene_places_no_label_on_a_node_at_both_window_canvas_sizes(name, canvas):
    width, height = canvas
    a = INPUTS[name]()
    scene, fig = _scene(a, render.compute_layout(a, width / height), (width / 100, height / 100))
    fig.canvas.draw()
    assert fig.get_size_inches() * fig.dpi == pytest.approx((width, height))
    assert elements(fig, "label")
    assert collisions(fig) == []
    # the plot fills the canvas: no disc lies off it
    assert all(inside(extent(fig, disc), fig) for disc in elements(fig, "node"))


LEFT_OUT = {  # input, canvas size: where nothing is cut, where both counts differ, where only some
    "none": ("lemon", SAVED_SIZE),
    "both": ("overlap", SMALL_SIZE),
    "only-left-out": ("three", SMALL_SIZE),
}


@pytest.mark.parametrize("case", LEFT_OUT)
def test_left_out_equals_the_counts_draw_figure_writes(case):
    name, size = LEFT_OUT[case]
    a = INPUTS[name]()
    layout = render.compute_layout(a, render.FIGURE_ASPECT)
    fresh, _ = drawn(a, size, layout=layout)
    written = [[int(n) for n in re.findall(r"\d+", t.get_text())] for t in frame(fresh, "left-out")]
    scene, _ = _scene_like(fresh, a, layout)
    if case == "none":
        assert written == []
        assert scene.left_out == (0, 0)
    else:
        (counts,) = written
        assert counts[0] != counts[1]  # shortened and left out cannot be swapped
        assert scene.left_out == tuple(counts)


def _discs(fig):
    """(letters, centre, radius) of every drawn node, in display pixels, measured on the drawing."""
    return [
        (
            disc.get_gid().removeprefix("node:"),
            fig.axes[0].transData.transform(disc.center),
            (extent(fig, disc)[2] - extent(fig, disc)[0]) / 2,
        )
        for disc in elements(fig, "node")
    ]


def _within(discs, point):
    """The node whose disc holds `point` (the nearest one): brute force."""
    found = [
        (math.dist(point, c), name) for name, c, radius in discs if math.dist(point, c) <= radius
    ]
    return min(found)[1] if found else None


def _node_at_defects(scene, points, expected):
    return [
        (p, node, scene.node_at(*p))
        for p, node in zip(points, expected)
        if scene.node_at(*p) != node
    ]


def _probe(name, resized=False):
    """A scene, points to ask about and the node under each, brute force on a drawn twin. With
    `resized` the canvas shrinks after the build and the scene is asked before any drawing."""
    a = SCENE_INPUTS[name]()
    width, height = CANVASES["default"]
    layout = render.compute_layout(a, width / height)

    def build():
        scene, fig = _scene(a, layout, (width / 100, height / 100))
        fig.canvas.draw()
        if resized:
            small_width, small_height = CANVASES["minimum"]
            fig.set_size_inches(small_width / 100, small_height / 100)
        return scene, fig

    _, twin = build()
    twin.canvas.draw()
    scene, _ = build()
    discs = _discs(twin)
    assert all(radius > 1 for _, _, radius in discs)
    centres = [centre for _, centre, _ in discs]
    between = [(p + q) / 2 for p, q in zip(centres, centres[1:])]
    inner = [c + (r - 1, 0) for _, c, r in discs]  # 1 px inside the drawn edge
    outer = [c + (r + 1, 0) for _, c, r in discs]  # 1 px outside it, on empty canvas
    points = [tuple(p) for p in (*centres, *between, *inner, *outer)] + [(1.0, 1.0)]
    expected = [_within(discs, p) for p in points]
    names = [name for name, _, _ in discs]
    assert expected[: len(discs)] == names
    assert expected[len(discs) * 2 - 1 : len(discs) * 3 - 1] == names  # inner points
    assert expected[len(discs) * 3 - 1 :] == [None] * (len(discs) + 1)  # outer points and a corner
    return scene, points, expected


PROBES = {"lemon": ("lemon", False), "overview": ("overview", False), "resized": ("lemon", True)}


@pytest.mark.parametrize("case", PROBES)
def test_node_at_returns_the_node_under_a_display_point(case):
    scene, points, expected = _probe(*PROBES[case])
    assert _node_at_defects(scene, points, expected) == []


def test_oracle_rejects_a_node_at_that_ignores_the_point():
    points, expected = _probe("lemon")[1:]

    class First:
        def node_at(self, x, y):
            return expected[0]

    assert _node_at_defects(First(), points, expected) != []


def test_hover_shows_the_letters_and_degrees():
    a = SCENE_INPUTS["overview"]()
    scene, fig = _scene(a, render.compute_layout(a, 1.5), (6, 4))
    (hover,) = elements(fig, "hover")
    assert not hover.get_visible()
    row = next(r for r in a.nodes if r.in_degree != r.out_degree)  # in and out cannot be swapped
    scene.show_hover(row.node)
    assert hover.get_visible() and row.node in hover.get_text()
    assert [int(n) for n in re.findall(r"\d+", hover.get_text())] == [row.in_degree, row.out_degree]
    scene.show_hover(None)
    assert not hover.get_visible()

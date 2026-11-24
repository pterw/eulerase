"""Tests for the layout and the figure: what is drawn, where, and what a saved file holds."""

import re

import pytest
from figure_probe import INPUTS, SAVED_SIZE, SMALL_SIZE, drawn, elements, extent, overlap
from genomes import random_dna
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.patches import Circle

from eulerase import render
from eulerase.analysis import Settings, analyze, settings_text
from eulerase.errors import EuleraseError
from eulerase.render import MAX_DRAW_NODES, compute_layout, draw_figure, save_figure
from eulerase.sequence import Record


def _analysis(sequence: str, k: int = 3, **settings):
    return analyze([Record("seq1", sequence)], Settings(k=k, **settings))


def _figure(a, upto=None):
    fig = Figure(figsize=(8, 6))
    layout = compute_layout(a)
    draw_figure(fig, a, layout, upto)
    return fig, layout


# A path with no repeated edge and no pair of opposite edges: ten distinct 3-mers in one trail.
PATH = "CCGTAATGCCTT"
MIN_DISC_PX = 8  # a node the eye can see


def _disc_defects(fig):
    """What is wrong with the node discs of a figure: one too small to see, or two that touch."""
    boxes = {disc.get_gid(): extent(fig, disc) for disc in elements(fig, "node")}
    found = [(gid, "small") for gid, box in boxes.items() if box[2] - box[0] < MIN_DISC_PX]
    gids = sorted(boxes)
    for i, first in enumerate(gids):
        found += [(first, g) for g in gids[i + 1 :] if overlap(boxes[first], boxes[g], margin=0)]
    return found


def test_oracle_rejects_invisible_and_touching_discs():
    fig = Figure(figsize=(4, 4))
    FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    for gid, centre, radius in (("node:AB", (20, 20), 5), ("node:BC", (60, 60), 0.5)):
        disc = Circle(centre, radius)
        disc.set_gid(gid)
        ax.add_patch(disc)
    fig.canvas.draw()
    assert _disc_defects(fig) == [("node:BC", "small")]
    touching = Circle((24, 20), 5)
    touching.set_gid("node:CD")
    ax.add_patch(touching)
    fig.canvas.draw()
    assert ("node:AB", "node:CD") in _disc_defects(fig)


@pytest.mark.parametrize("size", (SAVED_SIZE, SMALL_SIZE), ids=("saved", "small"))
@pytest.mark.parametrize("name", ("lemon", "repeat", "overlap"))
def test_every_node_is_drawn_visibly_and_apart(name, size):
    a = INPUTS[name]()
    fig, _ = drawn(a, size)
    assert len(elements(fig, "node")) == len(a.nodes)
    assert _disc_defects(fig) == []


def test_figure_title_carries_the_settings():
    settings = Settings(k=3, both_strands=True, min_multiplicity=1)
    a = analyze([Record("seq1", PATH)], settings)
    fig, _ = _figure(a)
    assert any(settings_text(settings) in t.get_text() for t in fig.texts)


def test_layout_is_none_above_the_drawing_limit_and_the_figure_says_why():
    a = _analysis(random_dna(1500, seed=1), k=7)
    assert len(a.nodes) > MAX_DRAW_NODES
    assert compute_layout(a) is None
    fig, _ = _figure(a)
    texts = [t.get_text() for t in elements(fig, "note")]
    assert len(texts) == len(fig.axes[0].texts) == 1
    assert str(MAX_DRAW_NODES) in texts[0] and str(len(a.nodes)) in texts[0]
    assert not elements(fig, "node")


# Four nodes, 78 arrows: the two counts of the drawing limit differ, so a swap shows.
REPEATED = "ACGT" * 20


def _drawing_counts(a):
    """The two counts of the drawing limit, read from the analysis: nodes, arrows."""
    return len(a.nodes), sum(len(t) for t in a.trails)


def test_oracle_rejects_a_node_count_used_as_the_arrow_count():
    assert _drawing_counts(_analysis(REPEATED)) == (4, 78)
    assert _drawing_counts(_analysis(REPEATED))[0] != 78


@pytest.mark.parametrize(
    "node_limit, arrow_limit, over",
    [
        (3, 91, {"nodes": 3}),
        (7, 50, {"arrows": 50}),
        (3, 50, {"nodes": 3, "arrows": 50}),
        (4, 78, {}),
    ],
    ids=["nodes only", "arrows only", "both", "exactly at both limits"],
)
def test_the_drawing_limit_counts_nodes_and_arrows(monkeypatch, node_limit, arrow_limit, over):
    monkeypatch.setattr(render, "MAX_DRAW_NODES", node_limit)
    monkeypatch.setattr(render, "MAX_DRAW_EDGES", arrow_limit)
    a = _analysis(REPEATED)
    nodes, arrows = _drawing_counts(a)
    assert nodes != arrows
    assert (nodes > node_limit, arrows > arrow_limit) == ("nodes" in over, "arrows" in over)

    fig, layout = _figure(a)
    if not over:
        assert layout is not None
        return
    assert layout is None
    (text,) = [t.get_text() for t in elements(fig, "note")]

    def mentions(number):
        return re.search(r"(?<!\d)" + str(number) + r"(?!\d)", text) is not None

    for name, count, limit in (("nodes", nodes, node_limit), ("arrows", arrows, arrow_limit)):
        if name in over:
            assert mentions(count) and mentions(limit)
        else:
            assert not mentions(limit)


FORMAT_STARTS = {
    ".png": (b"\x89PNG\r\n\x1a\n",),
    ".svg": (b"<?xml", b"<svg"),
    ".pdf": (b"%PDF",),
}


def test_save_figure_writes_each_format_and_refuses_another(tmp_path):
    a = _analysis(PATH)
    for suffix, starts in FORMAT_STARTS.items():
        path = tmp_path / f"graph{suffix}"
        save_figure(a, path)
        assert path.read_bytes().startswith(starts)
    written = sorted(p.name for p in tmp_path.iterdir())
    for refused in ("graph.bmp", "graph"):
        with pytest.raises(EuleraseError):
            save_figure(a, tmp_path / refused)
    assert sorted(p.name for p in tmp_path.iterdir()) == written

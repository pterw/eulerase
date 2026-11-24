"""Tests for the layout and the figure: what is drawn, where, and what a saved file holds."""

import math
import re

import pytest
from genomes import random_dna
from matplotlib.colors import to_rgb
from matplotlib.figure import Figure
from matplotlib.patches import Circle, FancyArrowPatch

from eulerase import render
from eulerase.analysis import Settings, analyze, settings_text
from eulerase.errors import EuleraseError
from eulerase.render import (
    EDGE,
    FUTURE,
    MAX_DRAW_NODES,
    compute_layout,
    draw_figure,
    save_figure,
)
from eulerase.sequence import Record

MIN_RADIUS, MAX_RADIUS = 0.01, 0.09
BRANCH_GROWTH = 1.3  # nodes that are branches, starts or ends are drawn larger above 40 nodes


def _analysis(sequence: str, k: int = 3, **settings):
    return analyze([Record("seq1", sequence)], Settings(k=k, **settings))


def _figure(a, upto=None):
    fig = Figure(figsize=(8, 6))
    layout = compute_layout(a)
    draw_figure(fig, a, layout, upto)
    return fig, layout


def _circles(fig):
    return [p for p in fig.axes[0].patches if isinstance(p, Circle)]


def _edge_labels(fig, a):
    kmers = {s.kmer for trail in a.trails for s in trail}
    return [t for t in fig.axes[0].texts if t.get_text().split(" ")[0] in kmers]


def _is_visible(radius: float, upper: float = MAX_RADIUS) -> bool:
    return MIN_RADIUS <= radius <= upper


# A path with no repeated edge and no pair of opposite edges: ten distinct 3-mers in one trail.
PATH = "CCGTAATGCCTT"


def test_oracle_rejects_zero_radius_nodes():
    assert _is_visible(0.05)
    assert not _is_visible(0.0)
    assert not _is_visible(0.001)


def test_every_node_is_drawn_visibly():
    a = _analysis(PATH)
    fig, _ = _figure(a)
    radii = [c.radius for c in _circles(fig)]
    assert len(radii) == len(a.nodes)
    assert all(_is_visible(r) for r in radii)


def test_small_graph_labels_every_edge_near_its_midpoint():
    a = _analysis(PATH)
    fig, layout = _figure(a)
    steps = [s for trail in a.trails for s in trail]
    labels = _edge_labels(fig, a)
    assert len(labels) == len(steps)
    midpoints = {
        s.kmer: tuple((layout[s.source][i] + layout[s.target][i]) / 2 for i in (0, 1))
        for s in steps
    }
    for label in labels:
        x, y = label.get_position()
        mx, my = midpoints[label.get_text().split(" ")[0]]
        assert math.hypot(x - mx, y - my) < 0.15


def test_dense_graph_keeps_branch_labels_and_drops_redundant_ones():
    a = _analysis(random_dna(120, seed=7), k=6)
    steps = [s for trail in a.trails for s in trail]
    assert len(steps) > 25 and len(a.nodes) > 40
    fig, layout = _figure(a)

    assert len(_edge_labels(fig, a)) < len(steps)

    preds, succs = {}, {}
    for s in steps:
        preds.setdefault(s.target, set()).add(s.source)
        succs.setdefault(s.source, set()).add(s.target)
    branches = {
        n.node for n in a.nodes if len(preds.get(n.node, ())) > 1 or len(succs.get(n.node, ())) > 1
    }
    assert branches
    labelled = {t.get_position() for t in fig.axes[0].texts}
    assert all(layout[n] in labelled for n in branches)
    assert len(labelled) < len(a.nodes) + len(steps)

    radii = [c.radius for c in _circles(fig)]
    assert len(radii) == len(a.nodes)
    assert all(_is_visible(r, MAX_RADIUS * BRANCH_GROWTH) for r in radii)


def test_a_step_marks_only_the_traversed_edges():
    a = _analysis(PATH)
    assert [len(t) for t in a.trails] == [10]
    fig, _ = _figure(a, upto=3)
    arrows = [p for p in fig.axes[0].patches if isinstance(p, FancyArrowPatch)]
    past = [p for p in arrows if p.get_edgecolor()[:3] == to_rgb(EDGE)]
    future = [p for p in arrows if p.get_edgecolor()[:3] == to_rgb(FUTURE)]
    assert (len(past), len(future)) == (3, 7)
    numbered = {t.get_text() for t in fig.axes[0].texts if "(#" in t.get_text()}
    first_three = {f"{s.kmer} (#{s.global_step})" for s in a.trails[0][:3]}
    assert numbered == first_three


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
    texts = [t.get_text() for t in fig.axes[0].texts]
    assert len(texts) == 1
    assert str(MAX_DRAW_NODES) in texts[0] and str(len(a.nodes)) in texts[0]
    assert not _circles(fig)


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
    (text,) = [t.get_text() for t in fig.axes[0].texts]

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

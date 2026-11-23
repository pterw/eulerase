"""The layout must not depend on the interpreter's string hashing (four processes, four seeds)."""

import itertools
import json
import math
import random
import subprocess
import sys

import pytest
from genomes import child_env

from eulerase.analysis import Settings, analyze
from eulerase.render import compute_layout
from eulerase.sequence import Record

CHILD = """
import json, random
import networkx as nx
from eulerase.analysis import Settings, analyze
from eulerase.render import compute_layout
from eulerase.sequence import Record

rng = random.Random(1)
first = "".join(rng.choice("ACGT") for _ in range(120))
second = "".join(rng.choice("LMNP") for _ in range(110))
a = analyze([Record("a", first), Record("b", second)], Settings(k=6))
g = nx.Graph()
g.add_nodes_from(n.node for n in a.nodes)
g.add_edges_from((s.source, s.target) for trail in a.trails for s in trail)
shape = (len(g), nx.number_connected_components(g))
assert shape == (209, 2), shape
print(json.dumps(sorted(compute_layout(a).items())))
"""


def _layouts_under(hash_seeds: list[int], scratch) -> list:
    """Start one child per hash seed together; return each child's layout, in seed order."""
    children = []
    try:
        for seed in hash_seeds:
            env = child_env()
            env["PYTHONHASHSEED"] = str(seed)
            env["MPLCONFIGDIR"] = str(scratch / f"mpl{seed}")
            children.append(
                subprocess.Popen(
                    [sys.executable, "-c", CHILD],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    env=env,
                )
            )
        layouts = []
        for child in children:
            out, err = child.communicate(timeout=120)
            assert child.returncode == 0, err
            layouts.append(json.loads(out))
        return layouts
    finally:
        for child in children:
            if child.poll() is None:
                child.kill()
                child.wait()


def _same_layout(one: list, other: list) -> bool:
    return one == other


def test_layout_spreads_two_components_inside_the_unit_square():
    rng = random.Random(1)
    first = "".join(rng.choice("ACGT") for _ in range(120))
    second = "".join(rng.choice("LMNP") for _ in range(110))
    a = analyze([Record("a", first), Record("b", second)], Settings(k=6))
    layout = compute_layout(a)
    # The two records share no letter, so a node's letters say which component it is in.
    boxes = []
    for alphabet in ("ACGT", "LMNP"):
        pts = [p for node, p in layout.items() if set(node) <= set(alphabet)]
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        boxes.append((min(xs), max(xs), min(ys), max(ys)))
    assert len(layout) == len(a.nodes) == 209
    assert len(set(layout.values())) == len(layout)
    # compute_layout centres the picture and scales its longer side to 2.0: [-1, 1] on each axis,
    # and the longer axis reaches both bounds (a layout scaled down would stay inside them).
    assert all(abs(c) <= 1.0 + 1e-9 for p in layout.values() for c in p)
    spans = [
        (min(p[i] for p in layout.values()), max(p[i] for p in layout.values())) for i in (0, 1)
    ]
    low, high = max(spans, key=lambda lo_hi: lo_hi[1] - lo_hi[0])
    assert low == pytest.approx(-1.0) and high == pytest.approx(1.0)
    # The layout follows the edges: an edge joins near points. Measured here: mean edge length
    # about 0.05 against a mean pair distance about 0.8 with the edges; about 0.5 against 0.9
    # when the layout ignores them. A quarter separates the two with room.
    ends = [(s.source, s.target) for trail in a.trails for s in trail if s.source != s.target]
    edge_mean = sum(math.dist(layout[u], layout[v]) for u, v in ends) / len(ends)
    pairs = list(itertools.combinations(layout.values(), 2))
    pair_mean = sum(math.dist(p, q) for p, q in pairs) / len(pairs)
    assert edge_mean < pair_mean / 4
    (ax0, ax1, ay0, ay1), (bx0, bx1, by0, by1) = boxes
    assert ax1 < bx0 or bx1 < ax0 or ay1 < by0 or by1 < ay0


def test_layout_is_identical_across_hash_seeds(tmp_path):
    first, *others = _layouts_under([1, 2, 3, 4], tmp_path)
    for other in others:
        assert _same_layout(first, other)


def test_layout_is_a_function_of_the_analysis_only():
    a = analyze([Record("a", "ACGTACGGTCA"), Record("b", "LMNPLMMNP")], Settings(k=4))
    assert compute_layout(a) == compute_layout(a)


def test_oracle_rejects_a_different_layout():
    assert _same_layout([["A", [0.0, 1.0]]], [["A", [0.0, 1.0]]])
    assert not _same_layout([["A", [0.0, 1.0]]], [["A", [0.0, 1.5]]])

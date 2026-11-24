"""The layout must not depend on the interpreter's string hashing (four processes, four seeds)."""

import dataclasses
import json
import math
import random
import subprocess
import sys

import pytest
from genomes import child_env

from eulerase.analysis import Settings, analyze
from eulerase.render import MAX_DRAW_NODES, compute_layout
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
print(json.dumps(sorted(compute_layout(a, 1.7).items())))
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


def _analysis(*sequences: str, k: int):
    return analyze([Record(f"r{i}", s) for i, s in enumerate(sequences)], Settings(k=k))


def _walk_order(a) -> list[str]:
    """The nodes in the order the walk first meets them, from the trails alone."""
    order: dict[str, None] = {}
    for trail in a.trails:
        order.setdefault(trail[0].source)
        for step in trail:
            order.setdefault(step.target)
    return list(order)


def _read_back(layout: dict) -> list[str]:
    """Nodes by row from the top, each row in its direction: even rows run left to right."""
    rows: dict[float, list[str]] = {}
    for node, (x, y) in layout.items():
        rows.setdefault(y, []).append(node)
    ordered = []
    for index, y in enumerate(sorted(rows, reverse=True)):
        ordered += sorted(rows[y], key=lambda n: layout[n][0], reverse=index % 2 == 1)
    return ordered


def _widest_columns(n: int, aspect: float) -> int:
    """Brute force: the column count with the widest spacing, the smaller on a tie."""
    best, best_spacing = 1, 0.0
    for c in range(1, n + 1):
        spacing = min(aspect / c, 1 / math.ceil(n / c))
        if spacing > best_spacing:
            best, best_spacing = c, spacing
    return best


def test_the_layout_matches_coordinates_written_in_the_test():
    # GTACG at k=3 walks GT, TA, AC, CG (not in sorted order).
    # Worked by hand from the rule: at aspect 2.0 the columns 2, 3 and 4 all give the spacing
    # 0.5, so the smaller, 2, is taken (a 2 by 2 grid whose second row runs right to left); at
    # aspect 4.0 four columns on one row give spacing 1.0.
    a = _analysis("GTACG", k=3)
    tie = {"GT": (0.5, 0.75), "TA": (1.5, 0.75), "AC": (1.5, 0.25), "CG": (0.5, 0.25)}
    one_row = {"GT": (0.5, 0.5), "TA": (1.5, 0.5), "AC": (2.5, 0.5), "CG": (3.5, 0.5)}
    assert len(set(tie.values())) == len(set(one_row.values())) == 4
    assert tie != one_row
    assert compute_layout(a, 2.0) == tie
    assert compute_layout(a, 4.0) == one_row


def _no_nodes():
    a = _analysis("ACGTA", k=3)
    return dataclasses.replace(a, nodes=(), trails=())


def _over_the_drawing_limit():
    rng = random.Random(5)
    a = _analysis("".join(rng.choice("ACGT") for _ in range(900)), k=6)
    assert len(a.nodes) > MAX_DRAW_NODES
    return a


def test_the_layout_has_no_nodes_for_no_input():
    assert compute_layout(_no_nodes(), 1.5) == {}


@pytest.mark.parametrize("aspect", [0.0, -1.0, math.nan, math.inf])
@pytest.mark.parametrize("case", ["a few nodes", "no nodes", "over the drawing limit"])
def test_a_bad_aspect_raises_before_anything_else_is_read(case, aspect):
    analysis = {
        "a few nodes": lambda: _analysis("ACGTA", k=3),
        "no nodes": _no_nodes,
        "over the drawing limit": _over_the_drawing_limit,
    }[case]()
    with pytest.raises(ValueError):
        compute_layout(analysis, aspect)


def test_the_grid_has_the_columns_with_the_widest_spacing():
    rng = random.Random(3)
    seen = set()
    for _ in range(60):
        text = "".join(rng.choice("ACGT") for _ in range(rng.randint(8, 70)))
        a = _analysis(text, k=rng.randint(3, 5))
        aspect = rng.uniform(0.5, 6.0)
        layout = compute_layout(a, aspect)
        top = max(y for _, y in layout.values())
        first_row = sum(1 for _, y in layout.values() if y == top)
        expected = _widest_columns(len(layout), aspect)
        assert first_row == expected
        seen.add(expected)
    assert len(seen) > 3  # the inputs reach several column counts, so a fixed count fails


@pytest.mark.parametrize(
    "sequences",
    [("ACGTTGCAACGTTGCA",), ("AACGGCAC", "GGACGCC", "AGAGGC")],
    ids=["one trail with a repeat", "three trails sharing nodes"],
)
def test_the_node_order_is_the_order_the_walk_first_meets_them(sequences):
    k = 4 if len(sequences) == 1 else 3
    a = _analysis(*sequences, k=k)
    order = _walk_order(a)
    assert len(a.trails) == len(sequences)
    assert order != sorted(order) and order != [n.node for n in a.nodes]
    assert len(order) == len(a.nodes) > 5
    assert _read_back(compute_layout(a, 1.0)) == order


@pytest.mark.parametrize("aspect", [0.5, 1.7, 6.0])
def test_layout_gives_every_node_one_distinct_point_inside_the_plot_area(aspect):
    rng = random.Random(1)
    first = "".join(rng.choice("ACGT") for _ in range(120))
    second = "".join(rng.choice("LMNP") for _ in range(110))
    a = _analysis(first, second, k=6)
    layout = compute_layout(a, aspect)
    assert len(layout) == len(a.nodes) == 209
    assert len(set(layout.values())) == len(layout)
    assert all(0 < x < aspect and 0 < y < 1 for x, y in layout.values())


def test_layout_is_identical_across_hash_seeds(tmp_path):
    first, *others = _layouts_under([1, 2, 3, 4], tmp_path)
    for other in others:
        assert _same_layout(first, other)


def test_layout_is_a_function_of_the_analysis_and_the_aspect_only():
    a = _analysis("ACGTACGGTCA", "LMNPLMMNP", k=4)
    again = _analysis("ACGTACGGTCA", "LMNPLMMNP", k=4)
    assert compute_layout(a, 1.5) == compute_layout(a, 1.5) == compute_layout(again, 1.5)
    # Rows of 2 at aspect 1.0 and one row of all the nodes at aspect 6.0.
    assert _widest_columns(len(a.nodes), 1.0) != _widest_columns(len(a.nodes), 6.0)
    assert compute_layout(a, 1.0) != compute_layout(a, 6.0)


def test_oracle_rejects_a_different_layout():
    assert _same_layout([["A", [0.0, 1.0]]], [["A", [0.0, 1.0]]])
    assert not _same_layout([["A", [0.0, 1.0]]], [["A", [0.0, 1.5]]])

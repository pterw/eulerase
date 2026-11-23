"""The layout must not depend on the interpreter's string hashing (two processes, two seeds)."""

import json
import subprocess
import sys

import pytest
from genomes import child_env

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


def _layout_under(hash_seed: int, scratch) -> list:
    env = child_env()
    env["PYTHONHASHSEED"] = str(hash_seed)
    env["MPLCONFIGDIR"] = str(scratch)
    done = subprocess.run(
        [sys.executable, "-c", CHILD],
        capture_output=True,
        text=True,
        timeout=120,
        env=env,
        check=True,
    )
    return json.loads(done.stdout)


def _same_layout(one: list, other: list) -> bool:
    return one == other


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="A4: node order follows string hashing",
)
def test_layout_is_identical_across_hash_seeds(tmp_path):
    one = _layout_under(1, tmp_path)
    two = _layout_under(2, tmp_path)
    assert _same_layout(one, two)


def test_oracle_rejects_a_different_layout():
    assert _same_layout([["A", [0.0, 1.0]]], [["A", [0.0, 1.0]]])
    assert not _same_layout([["A", [0.0, 1.0]]], [["A", [0.0, 1.5]]])

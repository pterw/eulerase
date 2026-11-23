"""Tests for the walk graph and the minimum Eulerian trail cover."""

import json
import random
import subprocess
import sys
from collections import Counter

import pytest
from genomes import brute_kmers, child_env

from eulerase import graph, kmers
from eulerase.errors import EuleraseError
from eulerase.sequence import Record


def _balances(items: list[str]) -> dict[str, int]:
    """Out minus in per (k-1)-mer, counted straight from the k-mers."""
    bal: Counter[str] = Counter()
    for km in items:
        bal[km[:-1]] += 1
        bal[km[1:]] -= 1
    return dict(bal)


def _minimum_trails(items: list[str]) -> int:
    """Spec 8.2: per weakly connected component the larger of 1 and the sum of the positive
    balances, summed over components. Components by union-find, independent of networkx."""
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for km in items:
        parent[find(km[:-1])] = find(km[1:])
    positive: Counter[str] = Counter()
    for node, b in _balances(items).items():
        positive[find(node)] += max(b, 0)
    return sum(max(1, positive[root]) for root in {find(n) for n in list(parent)})


def trail_errors(items: list[str], trails: list[graph.Trail]) -> list[str]:
    """Everything wrong with `trails` as a minimum cover of the walk graph of `items`."""
    errs = []
    used = Counter(edge for t in trails for edge in t)
    if any(c != 1 for c in used.values()):
        errs.append("an edge is used more than once")
    if Counter((u, v) for u, v, _ in used) != Counter((km[:-1], km[1:]) for km in items):
        errs.append("the edges used are not the edges of the input")
    for t in trails:
        if any(a[1] != b[0] for a, b in zip(t, t[1:])):
            errs.append("a trail is not a valid walk")
    k = len(items[0])
    spelled: Counter[str] = Counter()
    for t in trails:
        s = graph.spell(t)
        spelled.update(s[i : i + k] for i in range(len(s) - k + 1))
    if spelled != Counter(items):
        errs.append("the spelled k-mers differ from the input")
    if len(trails) != _minimum_trails(items):
        errs.append(f"{len(trails)} trails where {_minimum_trails(items)} are enough")
    return errs


def _random_items(rng: random.Random) -> list[str]:
    k = rng.randint(2, 6)
    alphabet = rng.choice(["AC", "ACG", "ACGT"])
    if rng.random() < 0.5:
        seqs = [
            "".join(rng.choice(alphabet) for _ in range(rng.randint(k, k + 12)))
            for _ in range(rng.randint(1, 4))
        ]
        both = alphabet == "ACGT" and rng.random() < 0.3
        items = brute_kmers(seqs, k, both)
    else:
        items = ["".join(rng.choice(alphabet) for _ in range(k)) for _ in range(rng.randint(1, 10))]
    return sorted(set(items)) if rng.random() < 0.4 else items


def _cover(items: list[str]) -> list[graph.Trail]:
    return graph.trail_cover(graph.build_walk_graph(items))


def test_trail_cover_properties_hold_on_random_inputs():
    rng = random.Random(2026)
    for _ in range(900):
        items = _random_items(rng)
        assert trail_errors(items, _cover(items)) == [], items


@pytest.mark.parametrize(
    "items",
    [
        ["AAA", "AAA", "AAC"],  # a self-loop with a parallel pair, then an exit
        ["ACG", "CGT", "CGA", "TTT"],  # two components, one needing two trails
    ],
)
def test_trail_cover_handles_self_loops_and_parallel_edges(items):
    assert trail_errors(items, _cover(items)) == []


_ORDER_PROGRAM = (
    "import json; from eulerase import graph; "
    "items = ['ACG', 'TCG', 'CGA', 'CCA', 'CCG', 'CCT', 'TCC', 'GGT', 'GTG', 'TGG']; "
    "cover = graph.trail_cover(graph.build_walk_graph(items)); "
    "print(json.dumps([graph.spell(t) for t in cover]))"
)


def test_trail_cover_order_is_the_same_in_every_process():
    # Two components: one with several sources and sinks, one balanced cycle.
    outputs = []
    for seed in ("1", "2"):
        env = child_env() | {"PYTHONHASHSEED": seed}
        run = subprocess.run(
            [sys.executable, "-c", _ORDER_PROGRAM],
            env=env,
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert run.returncode == 0, run.stderr
        outputs.append(json.loads(run.stdout))
    expected = ["ACG", "CCT", "TCCGA", "TCG", "CCA", "GGTGG"]
    assert outputs[0] == outputs[1] == expected


@pytest.mark.parametrize("fault", ["missing edge", "split trail"])
def test_oracle_rejects_a_missing_edge_and_a_split_trail(fault):
    items = ["ATG", "TGC", "GCA", "CAT"]
    (trail,) = _cover(items)
    assert trail_errors(items, [trail]) == []
    wrong = [trail[:-1]] if fault == "missing edge" else [trail[:2], trail[2:]]
    assert trail_errors(items, wrong) != []


def test_node_balance_matches_degree_difference():
    seqs = ["ATGCATGCA", "GGGAC", "TTTT"]
    items = brute_kmers(seqs, 3, both=True)
    expected = _balances(items)
    got = graph.node_balance(graph.build_walk_graph(items))
    assert got == {n: expected.get(n, 0) for n in got}
    assert set(got) == {km[:-1] for km in items} | {km[1:] for km in items}


def _items(sequence: str, k: int, *, times: int = 1, unique: bool = False) -> tuple[str, ...]:
    records = [Record(f"r{i}", sequence) for i in range(times)]
    return kmers.from_sequences(records, k, both_strands=False, unique=unique).items


def test_worked_case_a_read_given_ten_times():
    items = _items("ATGC", 3, times=10)
    g = graph.build_walk_graph(items)
    assert g.number_of_edges() == 20
    assert [graph.spell(t) for t in graph.trail_cover(g)] == ["ATGC"] * 10
    as_a_set = graph.trail_cover(graph.build_walk_graph(_items("ATGC", 3, times=10, unique=True)))
    assert [graph.spell(t) for t in as_a_set] == ["ATGC"]


def test_worked_case_a_repeat_inside_one_sequence():
    sequence = "ATGCATGCATGC"
    g = graph.build_walk_graph(_items(sequence, 3))
    assert g.number_of_edges() == 10
    assert Counter(d["kmer"] for _, _, d in g.edges(data=True)) == {
        "ATG": 3,
        "TGC": 3,
        "GCA": 2,
        "CAT": 2,
    }
    assert [graph.spell(t) for t in graph.trail_cover(g)] == [sequence]


@pytest.mark.parametrize("items", [[], ["ATG", "AT"], ["A", "C"]], ids=["empty", "mixed", "k=1"])
def test_walk_graph_refuses_bad_input(items):
    with pytest.raises(EuleraseError):
        graph.build_walk_graph(items)

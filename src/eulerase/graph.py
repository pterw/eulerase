"""The walk graph of a set of k-mers and its minimum cover by Eulerian trails."""

from collections import defaultdict
from collections.abc import Sequence

import networkx as nx

from eulerase.errors import EuleraseError

WalkEdge = tuple[str, str, int]
Trail = tuple[WalkEdge, ...]


def build_walk_graph(items: Sequence[str]) -> nx.MultiDiGraph:
    """Nodes = (k-1)-mers, edges = k-mers. Parallel edges for repeated k-mers."""
    if not items:
        raise EuleraseError("There are no k-mers to build a graph from; give some input.")
    lengths = {len(x) for x in items}
    if len(lengths) != 1:
        raise EuleraseError("All k-mers must have the same length; fix or remove the odd ones.")
    if lengths.pop() < 2:
        raise EuleraseError("k must be at least 2; give a larger k.")
    # Repeated k-mers become parallel edges. Copy counts live in the k-mer
    # multiset (see analysis), not on the edges.
    g = nx.MultiDiGraph()
    for km in items:
        g.add_edge(km[:-1], km[1:], kmer=km)
    return g


def node_balance(g: nx.MultiDiGraph) -> dict[str, int]:
    """Out-degree minus in-degree of every node."""
    return {n: g.out_degree(n) - g.in_degree(n) for n in g}


def _hierholzer(adj, start):
    """Iterative Hierholzer on adjacency lists of (v, key, is_dummy)."""
    stack = [(start, None)]
    out = []
    while stack:
        node, via = stack[-1]
        if adj[node]:
            v, key, dummy = adj[node].pop()
            stack.append((v, (node, v, key, dummy)))
        else:
            stack.pop()
            if via is not None:
                out.append(via)
    out.reverse()
    return out


def trail_cover(g: nx.MultiDiGraph) -> list[Trail]:
    """Minimum set of edge-disjoint trails covering every edge exactly once.
    One trail per component whenever an Eulerian path/circuit exists."""
    trails: list[Trail] = []
    for comp in sorted(nx.weakly_connected_components(g), key=min):
        bal = {n: g.out_degree(n) - g.in_degree(n) for n in comp}
        sources = [n for n in sorted(comp) for _ in range(max(bal[n], 0))]
        sinks = [n for n in sorted(comp) for _ in range(max(-bal[n], 0))]
        adj = defaultdict(list)
        for u, v, key in g.out_edges(comp, keys=True):
            adj[u].append((v, key, False))
        for i, (sk, sc) in enumerate(zip(sinks, sources)):
            adj[sk].append((sc, f"__dummy_{i}__", True))

        if not sources:
            start = min(comp, key=lambda n: (g.out_degree(n) == 0, n))
            circuit = _hierholzer(adj, start)
            real = tuple((u, v, k) for u, v, k, d in circuit if not d)
            if real:
                trails.append(real)
            continue

        raw = _hierholzer(adj, sinks[0])
        dummy_idx = [i for i, step in enumerate(raw) if step[3]]
        if not dummy_idx:
            trails.append(tuple((u, v, k) for u, v, k, d in raw if not d))
            continue
        first_d = dummy_idx[0]
        rotated = raw[first_d + 1 :] + raw[: first_d + 1]

        cur: list[WalkEdge] = []
        for u, v, key, dummy in rotated:
            if dummy:
                if cur:
                    trails.append(tuple(cur))
                    cur = []
            else:
                cur.append((u, v, key))
        if cur:
            trails.append(tuple(cur))
    return trails


def spell(trail: Trail) -> str:
    """Reconstruct sequence string spelled by an Eulerian trail."""
    if not trail:
        return ""
    return trail[0][0] + "".join(v[-1] for _, v, _ in trail)

"""Assembly graph, contigs and metrics: one edge per distinct k-mer, maximal paths as contigs."""

import statistics
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace

import networkx as nx

from eulerase.sequence import gc_content, reverse_complement


@dataclass(frozen=True)
class Contig:
    id: int
    sequence: str
    length: int
    edges: int
    kmer_count: int
    mean_multiplicity: float
    min_multiplicity: int
    gc: float | None
    circular: bool
    has_n: bool
    strand: str


@dataclass(frozen=True)
class Metrics:
    contigs: int
    total_length: int
    longest: int
    shortest: int
    mean_length: float | None
    n50: int
    l50: int
    n90: int
    gc: float | None
    mean_multiplicity: float | None


def build_assembly_graph(counts: Mapping[str, int]) -> nx.DiGraph:
    """One edge per distinct k-mer; `multiplicity` is how many times it was seen."""
    ag = nx.DiGraph()
    for km, n in counts.items():
        ag.add_edge(km[:-1], km[1:], kmer=km, multiplicity=n)
    return ag


def _least_rotation(s: str) -> str:
    """Lexicographically least rotation of *s* (Booth's algorithm, O(n))."""
    n = len(s)
    if n < 2:
        return s
    ss = s + s
    f = [-1] * (2 * n)
    least = 0
    for j in range(1, 2 * n):
        sj = ss[j]
        i = f[j - least - 1]
        while i != -1 and sj != ss[least + i + 1]:
            if sj < ss[least + i + 1]:
                least = j - i - 1
            i = f[i]
        if sj != ss[least + i + 1]:
            if sj < ss[least]:
                least = j
            f[j - least] = -1
        else:
            f[j - least] = i + 1
    return ss[least : least + n]


def _strand_key(contig: Contig) -> str:
    """Strand-independent identity of a contig: the smaller of its sequence and
    reverse complement. Circular contigs are compared up to rotation."""
    seq = contig.sequence
    if contig.circular:
        core = seq[: contig.edges]  # the cycle without its repeated k-1 tail
        return min(_least_rotation(core), _least_rotation(reverse_complement(core)))
    return min(seq, reverse_complement(seq))


def _drop_reverse_complement_twins(contigs: list[Contig]) -> list[Contig]:
    """Keep one contig per strand pair: the lexicographically smaller strand.
    A contig whose twin is absent (for example clipped) is kept as it is."""
    best: dict[tuple[bool, str], Contig] = {}
    for c in contigs:
        key = (c.circular, _strand_key(c))
        if key not in best or c.sequence < best[key].sequence:
            best[key] = c
    kept = {id(c) for c in best.values()}
    return [c for c in contigs if id(c) in kept]


def _make_contig(seq: str, multiplicities: list[int], circular: bool, strand: str) -> Contig:
    return Contig(
        id=0,
        sequence=seq,
        length=len(seq),
        edges=len(multiplicities),
        kmer_count=sum(multiplicities),
        mean_multiplicity=float(statistics.mean(multiplicities)),
        min_multiplicity=min(multiplicities),
        gc=gc_content(seq),
        circular=circular,
        has_n="N" in seq,
        strand=strand,
    )


def _open_contigs(ag: nx.DiGraph, strand: str, used_edges: set[tuple[str, str]]) -> list[Contig]:
    """Maximal non-branching paths from every branching or boundary node."""

    def is_internal(n: str) -> bool:
        return ag.in_degree(n) == 1 and ag.out_degree(n) == 1

    contigs: list[Contig] = []
    starts = sorted(n for n in ag if not is_internal(n) and ag.out_degree(n) > 0)
    for start in starts:
        for _, target in ag.out_edges(start):
            if (start, target) in used_edges:
                continue
            path = [(start, target)]
            used_edges.add((start, target))
            cur = target
            while is_internal(cur):
                succs = list(ag.successors(cur))
                if not succs or (cur, succs[0]) in used_edges:
                    break
                nxt = succs[0]
                path.append((cur, nxt))
                used_edges.add((cur, nxt))
                cur = nxt
            seq = path[0][0] + "".join(tgt[-1] for _, tgt in path)
            mults = [ag.edges[u, v]["multiplicity"] for u, v in path]
            contigs.append(_make_contig(seq, mults, False, strand))
    return contigs


def _cyclic_contigs(ag: nx.DiGraph, strand: str, used_edges: set[tuple[str, str]]) -> list[Contig]:
    """Isolated cycles: the edges no path used (every node on them is internal)."""
    contigs: list[Contig] = []
    remaining = set(ag.edges()) - used_edges
    if not remaining:
        return contigs
    temp = nx.DiGraph()
    for u, v in remaining:
        temp.add_edge(u, v, **ag.edges[u, v])
    for comp in nx.weakly_connected_components(temp):
        sub = temp.subgraph(comp)
        start_node = min(comp)
        cur, visited, cycle = start_node, set(), []
        while True:
            nxt_list = [w for w in sub.successors(cur) if (cur, w) not in visited]
            if not nxt_list:
                break
            nxt = nxt_list[0]
            cycle.append((cur, nxt))
            visited.add((cur, nxt))
            cur = nxt
            if cur == start_node:
                break
        if cycle:
            seq = cycle[0][0] + "".join(tgt[-1] for _, tgt in cycle)
            mults = [sub.edges[u, v]["multiplicity"] for u, v in cycle]
            contigs.append(_make_contig(seq, mults, True, strand))
    return contigs


def extract_contigs(ag: nx.DiGraph, k: int, *, both_strands: bool) -> tuple[Contig, ...]:
    """Contigs (maximal non-branching paths, then isolated cycles), longest first.

    With both strands the graph holds each strand, so each contig and its reverse
    complement count once: the lexicographically smaller strand is kept.

    Known limitation: the graph is collapsed, so repeats longer than k merge into one
    path; a pure repeat such as ``AAAAAAAAAA`` is reported as a short circular contig."""
    strand = "smaller-of-pair" if both_strands else "as-assembled"
    used_edges: set[tuple[str, str]] = set()
    contigs = _open_contigs(ag, strand, used_edges) + _cyclic_contigs(ag, strand, used_edges)
    if both_strands:
        contigs = _drop_reverse_complement_twins(contigs)
    contigs.sort(key=lambda c: -c.length)
    return tuple(replace(c, id=i) for i, c in enumerate(contigs, 1))


def _n_statistic(lengths: Sequence[int], threshold: float) -> tuple[int, int]:
    """Length and rank of the first contig, longest first, where the running total
    reaches *threshold*."""
    cumsum = 0
    for rank, length in enumerate(lengths, 1):
        cumsum += length
        if cumsum >= threshold:
            return length, rank
    raise ValueError("threshold above the total length")


def assembly_metrics(contigs: Sequence[Contig]) -> Metrics:
    """Standard assembly quality measures; with no contigs the integers are 0 and the rest None."""
    if not contigs:
        return Metrics(0, 0, 0, 0, None, 0, 0, 0, None, None)
    lengths = sorted((c.length for c in contigs), reverse=True)
    total = sum(lengths)
    n50, l50 = _n_statistic(lengths, total / 2)
    n90, _ = _n_statistic(lengths, total * 0.9)
    return Metrics(
        contigs=len(contigs),
        total_length=total,
        longest=lengths[0],
        shortest=lengths[-1],
        mean_length=total / len(contigs),
        n50=n50,
        l50=l50,
        n90=n90,
        gc=gc_content("".join(c.sequence for c in contigs)),
        mean_multiplicity=float(statistics.mean(c.mean_multiplicity for c in contigs)),
    )


def kmer_spectrum(counts: Mapping[str, int]) -> dict[int, int]:
    """Multiplicity to the number of distinct k-mers seen that many times, ascending."""
    return dict(sorted(Counter(counts.values()).items()))

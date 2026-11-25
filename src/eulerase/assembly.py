"""Assembly graph, contigs and metrics: one edge per distinct k-mer, maximal paths as contigs."""

import statistics
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from fractions import Fraction

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


def filter_min_multiplicity(ag: nx.DiGraph, minimum: int) -> tuple[nx.DiGraph, int]:
    """A copy without the k-mers seen fewer than *minimum* times, and how many distinct
    k-mers that removed. Nodes left without an edge are dropped too."""
    out = ag.copy()
    rare = [(u, v) for u, v, n in out.edges(data="multiplicity") if n < minimum]
    out.remove_edges_from(rare)
    out.remove_nodes_from(list(nx.isolates(out)))
    return out, len(rare)


def _chain(ag: nx.DiGraph, node: str, other: str, leaving: bool) -> tuple[tuple[str, str], ...]:
    """The branch that leaves *node* by the edge to *other* (or enters it by the edge from
    *other*): its edges in path order, ending at the first node that is not one edge in and
    one edge out."""
    step = ag.successors if leaving else ag.predecessors
    edges = []
    while True:
        edges.append((node, other) if leaving else (other, node))
        if ag.in_degree(other) != 1 or ag.out_degree(other) != 1:
            return tuple(edges if leaving else reversed(edges))
        node, other = other, next(iter(step(other)))


def clip_tips(ag: nx.DiGraph, k: int) -> tuple[nx.DiGraph, int, int]:
    """A copy without the dead-end branches that a stronger branch outweighs, the number of
    branches removed and the number of distinct k-mers (edges) that took.

    A branch is the chain of edges from a node to the next node that is not one edge in and
    one edge out. One that leaves its first node and ends where no edge leaves, or enters it
    and starts where no edge enters, is a dead end. It goes when it has fewer than ``2 * k``
    edges and another branch on the same side of the same first node has a strictly higher
    mean multiplicity. A chain is one branch however many of its ends make it a dead end.
    Every branch is judged on the graph as it comes in, and those that qualify go together."""
    doomed: set[tuple[tuple[str, str], ...]] = set()
    for node in ag:
        for leaving, degree, ends in (
            (True, ag.out_degree(node), ag.out_degree),
            (False, ag.in_degree(node), ag.in_degree),
        ):
            if degree < 2:
                continue
            others = ag.successors(node) if leaving else ag.predecessors(node)
            chains = [_chain(ag, node, other, leaving) for other in others]
            means = [Fraction(sum(ag.edges[e]["multiplicity"] for e in c), len(c)) for c in chains]
            best = max(means)
            for chain, mean in zip(chains, means, strict=True):
                far = chain[-1][1] if leaving else chain[0][0]
                if len(chain) < 2 * k and mean < best and ends(far) == 0:
                    doomed.add(chain)
    out = ag.copy()
    out.remove_edges_from(e for chain in doomed for e in chain)
    out.remove_nodes_from(list(nx.isolates(out)))
    return out, len(doomed), sum(len(c) for c in doomed)


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


def _opposite(contig: Contig) -> str:
    """The contig's reverse complement as it would be reported: a cycle is a circle, so the
    reverse complement is rotated to its least rotation, as the contig itself is."""
    opposite = reverse_complement(contig.sequence)
    return _least_rotation(opposite) if contig.circular else opposite


def _strand_key(contig: Contig) -> str:
    """Strand-independent identity of a contig: the smaller of its sequence and its
    reverse complement (a cycle's compared up to rotation)."""
    return min(contig.sequence, _opposite(contig))


def _is_own_twin(contig: Contig) -> bool:
    """True when the contig is its own reverse complement (up to rotation if a cycle)."""
    return contig.sequence == _opposite(contig)


def _drop_reverse_complement_twins(contigs: list[Contig]) -> list[Contig]:
    """Keep one contig per strand pair, the lexicographically smaller strand, labelled
    smaller-of-pair. A contig whose twin is absent (for example clipped) is kept as it is,
    labelled unpaired, unless it is its own reverse complement (then it is its own pair)."""
    groups: dict[tuple[bool, str], list[Contig]] = {}
    for c in contigs:
        groups.setdefault((c.circular, _strand_key(c)), []).append(c)
    kept = []
    for key, members in groups.items():
        best = min(members, key=lambda c: c.sequence)
        paired = len(members) > 1 or _is_own_twin(best)
        kept.append(replace(best, strand="smaller-of-pair" if paired else "unpaired"))
    return kept


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
                nxt = next(iter(ag.successors(cur)))
                path.append((cur, nxt))
                used_edges.add((cur, nxt))
                cur = nxt
            seq = path[0][0] + "".join(tgt[-1] for _, tgt in path)
            mults = [ag.edges[u, v]["multiplicity"] for u, v in path]
            contigs.append(_make_contig(seq, mults, False, strand))
    return contigs


def _cyclic_contigs(ag: nx.DiGraph, strand: str, used_edges: set[tuple[str, str]]) -> list[Contig]:
    """Isolated cycles: the edges no path used (every node on them is internal). A cycle of n
    k-mers is its circle written once: the first letter of each k-mer in cycle order, so n
    letters whatever k is. The walk starts at the smallest node, which makes the circle its
    least rotation (the smallest (k-1)-mer begins the smallest rotation)."""
    contigs: list[Contig] = []
    onward = {u: v for u, v in ag.edges() if (u, v) not in used_edges}  # one way on from each node
    for start in sorted(onward):
        if start not in onward:  # already walked as part of an earlier cycle
            continue
        cycle, node = [], start
        while node in onward:
            nxt = onward.pop(node)
            cycle.append((node, nxt))
            node = nxt
        seq = "".join(src[0] for src, _ in cycle)
        mults = [ag.edges[u, v]["multiplicity"] for u, v in cycle]
        contigs.append(_make_contig(seq, mults, True, strand))
    return contigs


def extract_contigs(ag: nx.DiGraph, k: int, *, both_strands: bool) -> tuple[Contig, ...]:
    """Contigs (maximal non-branching paths, then isolated cycles), longest first.

    With both strands the graph holds each strand, so each contig and its reverse
    complement count once: the lexicographically smaller strand is kept.

    Known limitation: the graph is collapsed, so repeats longer than k merge into one
    path; a pure repeat such as ``AAAAAAAAAA`` is reported as a cycle of one letter."""
    strand = "as-assembled"  # with both strands the twin pairing below relabels each contig
    used_edges: set[tuple[str, str]] = set()
    contigs = _open_contigs(ag, strand, used_edges) + _cyclic_contigs(ag, strand, used_edges)
    if both_strands:
        contigs = _drop_reverse_complement_twins(contigs)
    contigs.sort(key=lambda c: (-c.length, c.sequence))
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
        mean_multiplicity=sum(c.kmer_count for c in contigs) / sum(c.edges for c in contigs),
    )


def kmer_spectrum(counts: Mapping[str, int]) -> dict[int, int]:
    """Multiplicity to the number of distinct k-mers seen that many times, ascending."""
    return dict(sorted(Counter(counts.values()).items()))


Link = tuple[int, str, int, str]


def _flip(orient: str) -> str:
    return "-" if orient == "+" else "+"


def gfa_links(contigs: Sequence[Contig], k: int, both_strands: bool) -> list[Link]:
    """Links ``(from_id, from_orient, to_id, to_orient)`` between contigs.

    A link joins two oriented contigs that are not cycles when the last k-1 bases of the
    first equal the first k-1 of the second; a contig whose own ends overlap links to
    itself. With *both_strands* (each kept contig stands for its reverse complement too)
    either strand of every contig can take part, and ``-`` marks the reverse complement. A
    link and its mirror image (``b- -> a-`` for ``a+ -> b+``) are the same GFA edge, so only
    one is written. A cycle takes part in no overlap link: it has exactly one link, to
    itself, ``(id, "+", id, "+")``, with both strands too."""
    orients = ("+", "-") if both_strands else ("+",)
    seqs = {
        (c.id, o): c.sequence if o == "+" else reverse_complement(c.sequence)
        for c in contigs
        if not c.circular
        for o in orients
    }
    n = k - 1
    by_prefix: dict[str, list[tuple[int, str]]] = {}
    for key, seq in seqs.items():
        by_prefix.setdefault(seq[:n], []).append(key)
    links: list[Link] = []
    seen: set[Link] = set()
    for (a, ao), seq in seqs.items():
        for b, bo in by_prefix.get(seq[-n:], []):
            mirror = (b, _flip(bo), a, _flip(ao))
            if mirror in seen:
                continue
            link = (a, ao, b, bo)
            seen.add(link)
            links.append(link)
    return links + [(c.id, "+", c.id, "+") for c in contigs if c.circular]

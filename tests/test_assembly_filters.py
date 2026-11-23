"""Tests for the minimum-multiplicity filter, tip clipping, the contig partition and GFA links."""

import random
from collections import Counter

import pytest
from genomes import brute_kmers, random_dna, revcomp

from eulerase.assembly import (
    Contig,
    assembly_metrics,
    build_assembly_graph,
    clip_tips,
    extract_contigs,
    filter_min_multiplicity,
    gfa_links,
)
from eulerase.kmers import count_kmers


def _graph_kmers(ag) -> list[str]:
    return sorted(d["kmer"] for _, _, d in ag.edges(data=True))


def _contigs_of(reads: list[str], k: int, both: bool) -> tuple[Contig, ...]:
    counts = count_kmers(brute_kmers(reads, k, both))
    return extract_contigs(build_assembly_graph(counts), k, both_strands=both)


# --- minimum multiplicity -----------------------------------------------------------------


@pytest.mark.parametrize("minimum", [0, 1, 2, 3, 5, 1_000_000])
def test_min_multiplicity_removes_exactly_the_rarer_kmers(minimum):
    reads = [random_dna(60, seed=3), random_dna(60, seed=3), random_dna(40, seed=4)]
    items = brute_kmers(reads, 3, False)
    seen = Counter(items)
    ag = build_assembly_graph(count_kmers(items))
    kept, removed = filter_min_multiplicity(ag, minimum)
    expected = sorted(km for km, n in seen.items() if n >= minimum)
    assert removed == sum(1 for n in seen.values() if n < minimum)
    assert _graph_kmers(kept) == expected
    # No node is left behind without an edge, so a filter that removes everything gives an
    # empty graph.
    assert set(kept.nodes) == {x for km in expected for x in (km[:-1], km[1:])}


def test_removing_every_kmer_gives_an_empty_well_formed_result():
    ag = build_assembly_graph(count_kmers(brute_kmers([random_dna(30, seed=8)], 4, False)))
    kept, _ = filter_min_multiplicity(ag, 99)
    contigs = extract_contigs(kept, 4, both_strands=False)
    m = assembly_metrics(contigs)
    assert contigs == ()
    assert (m.contigs, m.total_length, m.n50) == (0, 0, 0)
    assert (m.mean_length, m.gc, m.mean_multiplicity) == (None, None, None)


@pytest.mark.parametrize(
    "apply",
    [lambda ag: filter_min_multiplicity(ag, 2), lambda ag: clip_tips(ag, 3)],
    ids=["filter_min_multiplicity", "clip_tips"],
)
def test_a_filter_returns_a_copy_and_leaves_its_input_alone(apply):
    counts = {"AAC": 3, "AAG": 1, "AGT": 1, "AGC": 2}
    ag = build_assembly_graph(counts)
    before = _graph_kmers(ag), sorted(ag.nodes)
    out = apply(ag)[0]
    assert out is not ag
    assert (_graph_kmers(ag), sorted(ag.nodes)) == before


# --- contigs partition the graph ------------------------------------------------------------


def _spelled(contigs: tuple[Contig, ...], k: int) -> list[tuple[int, str]]:
    """(contig index, k-mer) for every window a contig spells, from the sequences alone."""
    out = []
    for i, c in enumerate(contigs):
        windows = c.edges if c.circular else len(c.sequence) - k + 1
        out += [(i, c.sequence[j : j + k]) for j in range(windows)]
    return out


def _is_a_partition(contigs, ag, kept: set[str], k: int, both: bool) -> bool:
    """Every kept k-mer is spelled once (once per strand pair with both strands), and each
    contig is a path in the graph. `kept` comes from the counts, never from the contigs."""
    spelled = _spelled(contigs, k)
    names = [km for _, km in spelled]
    for c in contigs:
        walk = [c.sequence[j : j + k] for j in range(len(c.sequence) - k + 1)]
        if any(ag.get_edge_data(w[:-1], w[1:], {}).get("kmer") != w for w in walk):
            return False
        if c.circular and c.sequence[-(k - 1) :] != c.sequence[: k - 1]:
            return False
    if len(set(names)) != len(names) or not set(names) <= kept:
        return False
    if not both:
        return set(names) == kept
    # Both strands: a contig that is its own reverse complement spells a k-mer and its twin,
    # so the unit is the pair; each pair must be spelled, and by exactly one contig.
    owners: dict[frozenset[str], set[int]] = {}
    for i, km in spelled:
        owners.setdefault(frozenset((km, revcomp(km))), set()).add(i)
    pairs = {frozenset((km, revcomp(km))) for km in kept}
    return set(owners) == pairs and all(len(o) == 1 for o in owners.values())


def _random_case(case: int):
    rng = random.Random(case)
    alphabet = ("AC", "ACG", "ACGT")[case % 3]
    k = rng.randint(2, 6)
    both = bool(rng.getrandbits(1))
    reads = []
    for _ in range(rng.randint(1, 4)):
        chunk = "".join(rng.choice(alphabet) for _ in range(rng.randint(k, k + 12)))
        reads.append(chunk + chunk if rng.random() < 0.3 else chunk)
    if rng.random() < 0.5:
        reads.append(reads[0])
    return reads, k, both, rng.randint(1, 3)


@pytest.mark.parametrize("case", range(300))
def test_contigs_partition_the_assembly_graph_on_random_input(case):
    reads, k, both, minimum = _random_case(case)
    counts = count_kmers(brute_kmers(reads, k, both))
    kept = {km for km, n in counts.items() if n >= minimum}
    ag, _ = filter_min_multiplicity(build_assembly_graph(counts), minimum)
    contigs = extract_contigs(ag, k, both_strands=both)
    assert _is_a_partition(contigs, ag, kept, k, both)


def test_oracle_rejects_a_contig_set_that_repeats_a_kmer():
    k = 3
    counts = count_kmers(brute_kmers(["ACGTTGCA", "ACGTAA"], k, False))
    ag = build_assembly_graph(counts)
    contigs = extract_contigs(ag, k, both_strands=False)
    kept = set(counts)
    assert _is_a_partition(contigs, ag, kept, k, False)
    assert not _is_a_partition(contigs + contigs[:1], ag, kept, k, False)
    assert not _is_a_partition(contigs[:-1], ag, kept, k, False)


# --- hand-worked contigs: a branch, a self-loop, an isolated cycle -----------------------------


@pytest.mark.parametrize(
    ("reads", "k", "expected"),
    [
        # TACGG and TACGC at k=3: TAC and ACG twice, CGG and CGC once. Nodes TA, AC, CG, GG, GC;
        # edges TA>AC, AC>CG, CG>GG, CG>GC. CG has two edges out, so it ends the path from TA:
        # TACG (2 edges, 4 k-mers seen), then the two one-edge branches, CGC before CGG.
        (["TACGG", "TACGC"], 3, [(1, "TACG", 2, 4, False), (2, "CGC", 1, 1, False),
                                 (3, "CGG", 1, 1, False)]),
        # AAAA at k=3: AAA twice. One node AA with one edge to itself, so it is one cycle of
        # one edge: the sequence AA + A = AAA, circular, 2 k-mers seen.
        (["AAAA"], 3, [(1, "AAA", 1, 2, True)]),
        # ACGTAC at k=3: ACG, CGT, GTA, TAC once each, a ring AC>CG>GT>TA>AC of four edges; the
        # sequence is AC + G + T + A + C = ACGTAC (the first k-1 letters repeated at the end).
        (["ACGTAC"], 3, [(1, "ACGTAC", 4, 4, True)]),
    ],
    ids=["branch", "self-loop", "isolated-cycle"],
)
def test_contigs_of_hand_worked_graphs(reads, k, expected):
    got = _contigs_of(reads, k, False)
    assert [(c.id, c.sequence, c.edges, c.kmer_count, c.circular) for c in got] == expected


# --- tip clipping ---------------------------------------------------------------------------

_A2 = pytest.mark.xfail(
    strict=True, raises=AssertionError, reason="A2: tip clipping compares nothing"
)

_LONG = {  # two branches of 6 = 2k edges from the node AA at k=3; the first is seen once
    **{km: 1 for km in ("AAC", "ACG", "CGT", "GTC", "TCA", "CAT")},
    **{km: 5 for km in ("AAG", "AGG", "GGC", "GCC", "CCT", "CTT")},
}


@pytest.mark.parametrize(
    ("counts", "survivors", "branches", "kmers"),
    [
        # u=AA, p=AC, b=AG, x=GT, y=GC. Leaving u: p (mean 3) and b (mean 3) tie, both kept.
        # Leaving b: x (1) against y (2): only x goes.
        pytest.param(
            {"AAC": 3, "AAG": 3, "AGT": 1, "AGC": 2},
            ["AAC", "AAG", "AGC"],
            1,
            1,
            id="loses-only-AGT",
            marks=_A2,
        ),
        pytest.param(
            {"AAC": 10, "AAG": 1}, ["AAC"], 1, 1, id="seen-10-times-beats-seen-once", marks=_A2
        ),
        pytest.param(
            {"AAC": 2, "AAG": 2}, ["AAC", "AAG"], 0, 0, id="equal-branches-both-kept", marks=_A2
        ),
        pytest.param(_LONG, sorted(_LONG), 0, 0, id="branch-of-2k-edges-is-kept"),
        # ACAAAATA at k=3: ACA, CAA, AAA twice, AAT, ATA. The self-loop AAA (mean 2) beats both
        # the tail AAT, ATA (mean 1) and the head ACA, CAA (mean 1): four k-mers, two branches.
        pytest.param(
            {"ACA": 1, "CAA": 1, "AAA": 2, "AAT": 1, "ATA": 1}, ["AAA"], 2, 4, id="self-loop-wins"
        ),
    ],
)
def test_tip_clipping_follows_the_rule(counts, survivors, branches, kmers):
    clipped, n_branches, n_kmers = clip_tips(build_assembly_graph(counts), 3)
    assert (_graph_kmers(clipped), n_branches, n_kmers) == (survivors, branches, kmers)


@_A2
def test_tip_clipping_gives_the_same_contigs_whatever_order_the_reads_came_in():
    # The mapping is built in read order, so the graph's node order differs between the two.
    # Each way leaves one contig, so the contig numbering is not in question.
    a, b = ["TACGT"] * 4, ["TACGA"]
    outcomes = []
    for reads in (a + b, b + a):
        counts = dict(Counter(brute_kmers(reads, 3, False)))
        clipped, _, _ = clip_tips(build_assembly_graph(counts), 3)
        outcomes.append(extract_contigs(clipped, 3, both_strands=False))
    assert outcomes[0] == outcomes[1]


# --- GFA links ------------------------------------------------------------------------------


def _oriented(contigs, cid: int, orient: str) -> str:
    seq = next(c.sequence for c in contigs if c.id == cid)
    return seq if orient == "+" else revcomp(seq)


def _branching_reads(seed: int) -> list[str]:
    stem, tail1, tail2, head = (random_dna(40, seed * 3 + i) for i in (0, 1, 2, 100))
    return [head + stem + tail1, head + stem + tail2]


@pytest.mark.parametrize("seed", range(3))
def test_gfa_links_branching_both_strands(seed):
    k = 11
    contigs = _contigs_of(_branching_reads(seed), k, True)
    links = gfa_links(contigs, k, True)
    assert len(contigs) == 3  # the shared head and stem, and the two tails, once each
    assert len(links) == 2
    for a, ao, b, bo in links:
        assert _oriented(contigs, a, ao)[-(k - 1) :] == _oriented(contigs, b, bo)[: k - 1]


def test_gfa_links_forward_only_are_plus():
    k = 11
    stem, tail1, tail2 = (random_dna(40, 60 + i) for i in range(3))
    contigs = _contigs_of([stem + tail1, stem + tail2], k, False)
    links = gfa_links(contigs, k, False)
    assert len(links) == 2
    for a, ao, b, bo in links:
        assert (ao, bo) == ("+", "+")
        assert _oriented(contigs, a, ao)[-(k - 1) :] == _oriented(contigs, b, bo)[: k - 1]

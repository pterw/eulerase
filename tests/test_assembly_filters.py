"""Tests for the minimum-multiplicity filter, tip clipping, the contig partition and GFA links."""

import random
from collections import Counter
from fractions import Fraction

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
    """(contig index, k-mer) for every window a contig spells, from the sequences alone: a line
    has len - k + 1 windows, a cycle (its circle written once) has n, read round the circle."""
    out = []
    for i, c in enumerate(contigs):
        if c.circular:
            ring = c.sequence * (k // len(c.sequence) + 2)
            out += [(i, ring[j : j + k]) for j in range(len(c.sequence))]
        else:
            out += [(i, c.sequence[j : j + k]) for j in range(len(c.sequence) - k + 1)]
    return out


def _is_a_partition(contigs, ag, kept: set[str], k: int, both: bool) -> bool:
    """Every kept k-mer is spelled once (once per strand pair with both strands), and each
    contig is a path in the graph. `kept` comes from the counts, never from the contigs."""
    spelled = _spelled(contigs, k)
    names = [km for _, km in spelled]
    if any(ag.get_edge_data(w[:-1], w[1:], {}).get("kmer") != w for w in names):
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
        # one edge: the circle is its one letter A, circular, 2 k-mers seen.
        (["AAAA"], 3, [(1, "A", 1, 2, True)]),
        # ACGTAC at k=3: ACG, CGT, GTA, TAC once each, a ring AC>CG>GT>TA>AC of four edges; the
        # circle is the first letter of each k-mer, A C G T, which is already its least rotation.
        (["ACGTAC"], 3, [(1, "ACGT", 4, 4, True)]),
    ],
    ids=["branch", "self-loop", "isolated-cycle"],
)
def test_contigs_of_hand_worked_graphs(reads, k, expected):
    got = _contigs_of(reads, k, False)
    assert [(c.id, c.sequence, c.edges, c.kmer_count, c.circular) for c in got] == expected


# --- tip clipping ---------------------------------------------------------------------------

_LONG = {  # two branches of 6 = 2k edges from the node AA at k=3; the first is seen once
    **{km: 1 for km in ("AAC", "ACG", "CGT", "GTC", "TCA", "CAT")},
    **{km: 5 for km in ("AAG", "AGG", "GGC", "GCC", "CCT", "CTT")},
}
# The same two branches with the letters of every k-mer reversed, so that they enter AA: the
# weak one runs TA>AC>CT>TG>GC>CA>AA (TAC, ACT, CTG, TGC, GCA, CAA), the strong one the like.
_LONG_IN = {km[::-1]: n for km, n in _LONG.items()}
# The weak branch cut to 5 = 2k - 1 edges (the last k-mer CAT dropped), the strong one left at 6.
_SHORT = {km: n for km, n in _LONG.items() if km != "CAT"}
_SHORT_IN = {km[::-1]: n for km, n in _SHORT.items()}


def _strong(counts):
    return sorted(km for km, n in counts.items() if n == 5)


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
        ),
        pytest.param({"AAC": 10, "AAG": 1}, ["AAC"], 1, 1, id="seen-10-times-beats-seen-once"),
        pytest.param({"AAC": 2, "AAG": 2}, ["AAC", "AAG"], 0, 0, id="equal-branches-both-kept"),
        # Leaving AA: AAC and AAG (once each) tie, AAT (5) is strictly better than both: both go.
        pytest.param(
            {"AAC": 1, "AAG": 1, "AAT": 5}, ["AAT"], 2, 2, id="two-equal-branches-beside-a-better"
        ),
        # Leaving AA: AAC>ACG (2 edges, seen 2 + 2 = 4 times, mean 2) against AAT (1 edge, seen
        # 3 times, mean 3). AAT has the smaller total and the higher mean: the two-edge branch goes.
        pytest.param(
            {"AAC": 2, "ACG": 2, "AAT": 3}, ["AAT"], 1, 2, id="the-mean-decides-not-the-total"
        ),
        pytest.param(_LONG, sorted(_LONG), 0, 0, id="branch-of-2k-edges-is-kept-leaving"),
        pytest.param(_LONG_IN, sorted(_LONG_IN), 0, 0, id="branch-of-2k-edges-is-kept-entering"),
        # 5 = 2k - 1 edges is short enough: the weak branch goes, the 6 strong edges stay.
        pytest.param(_SHORT, _strong(_SHORT), 1, 5, id="branch-of-2k-minus-1-edges-goes-leaving"),
        pytest.param(
            _SHORT_IN, _strong(_SHORT_IN), 1, 5, id="branch-of-2k-minus-1-edges-goes-entering"
        ),
        # ACAAAATA at k=3: ACA, CAA, AAA twice, AAT, ATA. The self-loop AAA (mean 2) beats both
        # the tail AAT, ATA (mean 1) and the head ACA, CAA (mean 1): four k-mers, two branches.
        pytest.param(
            {"ACA": 1, "CAA": 1, "AAA": 2, "AAT": 1, "ATA": 1}, ["AAA"], 2, 4, id="self-loop-wins"
        ),
        # AAC runs AA>AC. AA has no edge in and two out; AC has two in and none out; TA has none
        # in. Leaving AA the chain AAC is a dead end with AAG (5) beside it; entering AC it is
        # a dead end with TAC (5) beside it. One chain, so one branch and one k-mer.
        pytest.param(
            {"AAC": 1, "AAG": 5, "TAC": 5},
            ["AAG", "TAC"],
            1,
            1,
            id="a-chain-dead-at-both-ends-is-one-branch",
        ),
    ],
)
def test_tip_clipping_follows_the_rule(counts, survivors, branches, kmers):
    clipped, n_branches, n_kmers = clip_tips(build_assembly_graph(counts), 3)
    assert (_graph_kmers(clipped), n_branches, n_kmers) == (survivors, branches, kmers)


def _random_counts(seed: int) -> tuple[dict[str, int], int]:
    """k-mer counts of 1 to 4 random reads over A and C, each read repeated 1 to 3 times."""
    rng = random.Random(seed)
    k = rng.randint(3, 5)
    reads = []
    for _ in range(rng.randint(1, 4)):
        read = "".join(rng.choice("AC") for _ in range(rng.randint(k, k + 10)))
        reads += [read] * rng.randint(1, 3)
    return dict(count_kmers(brute_kmers(reads, k, False))), k


def _judge(counts: dict[str, int], k: int, *, by_mean: bool = True) -> dict:
    """The rule's sentences applied literally to the k-mer counts, with no graph library.

    A branch is a chain of k-mers from a first node. It leaves the node (its first k-mer
    starts there) or enters it (its last k-mer ends there), and runs through nodes with one
    k-mer in and one out up to the next node that is not like that. A chain is identified by
    its ordered k-mers, so one judged from both of its ends is one branch. Returns the removed
    chains, the chains that are a short dead end from both ends, and how many short dead ends
    have an equal rival and were kept. `by_mean=False` is a wrong judge for which any other
    branch on the same side counts as better."""
    into: dict[str, list[str]] = {}
    out_of: dict[str, list[str]] = {}
    for km in counts:
        out_of.setdefault(km[:-1], []).append(km)
        into.setdefault(km[1:], []).append(km)

    def grow(first: str, leaving: bool) -> tuple[str, ...]:
        chain = [first]
        while True:
            node = chain[-1][1:] if leaving else chain[0][:-1]
            if len(into.get(node, [])) != 1 or len(out_of.get(node, [])) != 1:
                return tuple(chain)
            step = out_of[node][0] if leaving else into[node][0]
            if step in chain:  # a lone cycle
                return tuple(chain)
            chain = [*chain, step] if leaving else [step, *chain]

    branches = []  # (first node, leaving?, chain, mean, a short dead end?)
    for node in sorted(set(into) | set(out_of)):
        for leaving in (True, False):
            for first in (out_of if leaving else into).get(node, []):
                chain = grow(first, leaving)
                far = chain[-1][1:] if leaving else chain[0][:-1]
                dead = not (out_of if leaving else into).get(far, [])
                mean = Fraction(sum(counts[km] for km in chain), len(chain))
                branches.append((node, leaving, chain, mean, dead and len(chain) < 2 * k))

    def rivals(b):
        return [o for o in branches if o is not b and o[:2] == b[:2]]

    def goes(b) -> bool:
        return b[4] and any((o[3] > b[3]) if by_mean else True for o in rivals(b))

    short_dead = [b[2] for b in branches if b[4]]
    return {
        "removed": {b[2] for b in branches if goes(b)},
        "both_ends": {c for c in set(short_dead) if short_dead.count(c) == 2},
        "tied_and_kept": sum(
            1 for b in branches if b[4] and not goes(b) and any(o[3] == b[3] for o in rivals(b))
        ),
    }


def test_tip_clipping_matches_a_brute_force_judge():
    saw_clipped = saw_tie_kept = saw_both_ends = saw_counts_differ = False
    for seed in range(400):
        counts, k = _random_counts(seed)
        verdict = _judge(counts, k)
        clipped, n_branches, n_kmers = clip_tips(build_assembly_graph(counts), k)
        removed = {km for chain in verdict["removed"] for km in chain}
        assert _graph_kmers(clipped) == sorted(set(counts) - removed), (seed, counts, k)
        assert (n_branches, n_kmers) == (len(verdict["removed"]), len(removed)), (seed, counts)
        saw_clipped |= bool(removed)
        saw_tie_kept |= verdict["tied_and_kept"] > 0
        saw_both_ends |= bool(verdict["both_ends"])
        saw_counts_differ |= n_branches != n_kmers  # else branches and k-mers could be swapped
    assert (saw_clipped, saw_tie_kept, saw_both_ends, saw_counts_differ) == (True,) * 4


def test_oracle_rejects_a_judge_that_ignores_multiplicity():
    # AAC seen 10 times, AAG once, both dead ends of AA: only AAG is worse, so only AAG goes.
    counts = {"AAC": 10, "AAG": 1}
    assert counts["AAC"] != counts["AAG"]  # the rivals must differ or the two judges agree
    assert _judge(counts, 3)["removed"] == {("AAG",)}
    assert _judge(counts, 3, by_mean=False)["removed"] == {("AAC",), ("AAG",)}


def test_clipping_is_independent_of_read_order():
    # The graph's node and edge order follows the order of the mapping, so shuffling it
    # changes the order in which the clipper meets every node.
    saw_clipped = False
    for seed in range(60):
        counts, k = _random_counts(seed)
        items = list(counts.items())
        outcomes = []
        for shuffle_seed in range(4):
            random.Random(shuffle_seed).shuffle(items)
            clipped, branches, kmers = clip_tips(build_assembly_graph(dict(items)), k)
            contigs = extract_contigs(clipped, k, both_strands=False)
            outcomes.append((_graph_kmers(clipped), branches, kmers, contigs))
        assert all(o == outcomes[0] for o in outcomes), (seed, counts)
        saw_clipped |= outcomes[0][1] > 0
    assert saw_clipped


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

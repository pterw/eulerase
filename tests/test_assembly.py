"""Tests for the assembly graph, contigs, metrics and the k-mer spectrum."""

import json
import os
import random
import statistics
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest
from genomes import brute_kmers, child_env, is_the_genome, least_rotation, random_dna, revcomp
from pipeline import contigs_of

from eulerase.assembly import (
    assembly_metrics,
    build_assembly_graph,
    extract_contigs,
    filter_min_multiplicity,
    gfa_links,
    kmer_spectrum,
)
from eulerase.kmers import count_kmers


def _repeat_free_genome(length: int, k: int, both: bool, seed: int) -> str:
    """A genome with no repeated (k-1)-mer and, with both strands, none shared between them."""
    n = k - 1
    for attempt in range(2000):
        g = random_dna(length, seed * 10_000 + attempt)
        words = [g[i : i + n] for i in range(len(g) - n + 1)]
        if both:
            r = revcomp(g)
            words += [r[i : i + n] for i in range(len(r) - n + 1)]
        if len(set(words)) == len(words):
            return g
    raise AssertionError("no repeat-free genome found")


def _tiles(genome: str, size: int = 80, step: int = 40) -> list[str]:
    """Overlapping reads that cover the genome; neighbours share at least 40 letters."""
    return [genome[i : i + size] for i in range(0, len(genome) - size + 1, step)] + [genome[-size:]]


@pytest.mark.parametrize("both", [False, True])
@pytest.mark.parametrize("seed", range(20))
def test_error_free_genome_gives_one_contig_equal_to_it(seed, both):
    k = 11 + seed % 15
    genome = _repeat_free_genome(300 + (seed * 37) % 701, k, both, seed)
    contigs = contigs_of(_tiles(genome), k, both)
    assert is_the_genome([c.sequence for c in contigs], genome, both)
    assert contigs[0].strand == ("smaller-of-pair" if both else "as-assembled")


@pytest.mark.parametrize(
    ("genome", "k"),
    [
        (random_dna(3000, seed=11), 21),
        (random_dna(400, seed=5), 15),
        ("CTTCGCCTGATACG", 5),
    ],
)
def test_both_strands_gives_one_canonical_contig_per_strand_pair(genome, k):
    contigs = contigs_of([genome], k, True)
    seqs = [c.sequence for c in contigs]
    assert len(genome) - 5 <= sum(len(s) for s in seqs) <= len(genome) + 5
    for s in seqs:
        assert revcomp(s) not in seqs or revcomp(s) == s
        assert s == min(s, revcomp(s))


def test_forward_input_is_never_deduplicated():
    read = "CTTCGCCTGATACG"
    with_twin = contigs_of([read, revcomp(read)], 5, False)
    assert sorted(c.sequence for c in with_twin) == sorted([read, revcomp(read)])
    unrelated = [random_dna(40, seed=21), random_dna(40, seed=22)]
    assert sorted(c.sequence for c in contigs_of(unrelated, 9, False)) == sorted(unrelated)


def test_n_does_not_defeat_twin_removal():
    read = "ACGTNACGTTGCAGGATCCATTGACN"
    both = contigs_of([read], 5, True)
    seqs = [c.sequence for c in both]
    assert all(revcomp(s) not in seqs or revcomp(s) == s for s in seqs)
    assert any("N" in s for s in seqs)
    both_strand_graph = build_assembly_graph(count_kmers(brute_kmers([read], 5, True)))
    undeduplicated = extract_contigs(both_strand_graph, 5, both_strands=False)
    assert len(both) < len(undeduplicated)
    # has_n follows the letters of each contig, in both directions.
    mixed = contigs_of([read, random_dna(30, seed=31)], 7, False)
    assert {c.has_n for c in mixed} == {True, False}
    assert all(c.has_n == ("N" in c.sequence) for c in mixed)


def test_circular_genome_with_both_strands_is_one_circular_contig():
    genome = random_dna(300, seed=9)
    k = 15
    contigs = contigs_of([genome + genome[: k - 1]], k, True)
    assert len(contigs) == 1
    assert contigs[0].circular is True
    assert contigs[0].edges == len(genome)
    assert contigs[0].length == len(genome)
    assert contigs[0].sequence in {least_rotation(genome), least_rotation(revcomp(genome))}


def test_multiplicity_is_per_distinct_kmer():
    ag = build_assembly_graph({"AAA": 3, "AAC": 1})
    assert ag.edges["AA", "AA"]["multiplicity"] == 3
    assert ag.edges["AA", "AC"]["multiplicity"] == 1


def _n_statistic(lengths: list[int], fraction: float) -> tuple[int, int]:
    """Independent N50-style answer: (length, count) of the longest contigs that
    first hold at least `fraction` of the total."""
    ordered = sorted(lengths, reverse=True)
    for count in range(1, len(ordered) + 1):
        if sum(ordered[:count]) >= fraction * sum(ordered):
            return ordered[count - 1], count
    raise AssertionError("unreachable")


def _windows(sequence: str, k: int, circular: bool) -> list[str]:
    """The k-mers a contig spells: the windows of a line, or the n of a circle read round."""
    if circular:
        ring = sequence * (k // len(sequence) + 2)
        return [ring[j : j + k] for j in range(len(sequence))]
    return [sequence[j : j + k] for j in range(len(sequence) - k + 1)]


def test_oracle_rejects_a_wrong_n50():
    assert _n_statistic([10, 10, 10, 1], 0.5) == (10, 2)
    assert _n_statistic([10, 10, 10, 1], 0.5) != (1, 4)


@pytest.mark.parametrize("seed", range(150))
def test_metrics_match_direct_computation(seed):
    k = 4 + seed % 12
    reads = [random_dna(k + 5 + (seed * (i + 3)) % 40, seed * 7 + i) for i in range(1 + seed % 5)]
    counts = count_kmers(brute_kmers(reads, k, False) + brute_kmers(reads[:1], k, False))
    contigs = extract_contigs(build_assembly_graph(counts), k, both_strands=False)
    windows = total_seen = 0
    for c in contigs:
        seen = [counts[w] for w in _windows(c.sequence, k, c.circular)]
        assert (c.edges, c.kmer_count, c.min_multiplicity) == (len(seen), sum(seen), min(seen))
        assert c.mean_multiplicity == pytest.approx(sum(seen) / len(seen))
        assert c.gc == (c.sequence.count("G") + c.sequence.count("C")) / len(c.sequence)
        windows += len(seen)
        total_seen += sum(seen)
    lengths = [len(c.sequence) for c in contigs]
    joined = "".join(c.sequence for c in contigs)
    m = assembly_metrics(contigs)
    assert (m.contigs, m.total_length) == (len(contigs), sum(lengths))
    assert (m.longest, m.shortest) == (max(lengths), min(lengths))
    assert m.mean_length == pytest.approx(sum(lengths) / len(lengths))
    assert (m.n50, m.l50) == _n_statistic(lengths, 0.5)
    assert m.n90 == _n_statistic(lengths, 0.9)[0]
    assert m.gc == pytest.approx((joined.count("G") + joined.count("C")) / len(joined))
    assert m.mean_multiplicity == pytest.approx(total_seen / windows)


def test_no_contigs_gives_zero_integers_and_missing_values():
    contigs = extract_contigs(build_assembly_graph({}), 3, both_strands=False)
    assert contigs == ()
    m = assembly_metrics(contigs)
    assert (m.contigs, m.total_length, m.longest, m.shortest, m.n50, m.l50, m.n90) == (0,) * 7
    assert (m.mean_length, m.gc, m.mean_multiplicity) == (None, None, None)


def test_kmer_spectrum_matches_counter_of_counts():
    for seed in range(100):
        items = [random_dna(3, seed * 31 + i) for i in range(5 + seed % 40)]
        expected = {
            m: n
            for m in range(1, len(items) + 1)
            if (n := sum(1 for x in set(items) if items.count(x) == m))
        }
        spectrum = kmer_spectrum(count_kmers(items))
        assert list(spectrum.items()) == list(expected.items())


def test_contigs_are_numbered_in_a_stable_order():
    # Reads ACGTTGCACCATGGATCC and ACGTTGCAGGTACCTAGG at k=3 hold 24 distinct 3-mers.
    # The 2-mers CG, TT, GC, TC, CT and GA each have one edge in and one out, so they join
    # their neighbours: ACG+CGT = ACGT, GTT+TTG = GTTG, TGC+GCA = TGCA, ATC+TCC = ATCC,
    # CCT+CTA = CCTA, GGA+GAT = GGAT (six contigs of length 4). Every other 2-mer branches,
    # so the remaining 12 3-mers are contigs of length 3. Longest first, ties by sequence,
    # ids from 1, whatever the read order.
    expected = [
        (1, 4, "ACGT"),
        (2, 4, "ATCC"),
        (3, 4, "CCTA"),
        (4, 4, "GGAT"),
        (5, 4, "GTTG"),
        (6, 4, "TGCA"),
        (7, 3, "ACC"),
        (8, 3, "AGG"),
        (9, 3, "ATG"),
        (10, 3, "CAC"),
        (11, 3, "CAG"),
        (12, 3, "CAT"),
        (13, 3, "CCA"),
        (14, 3, "GGT"),
        (15, 3, "GTA"),
        (16, 3, "TAC"),
        (17, 3, "TAG"),
        (18, 3, "TGG"),
    ]
    reads = ["ACGTTGCACCATGGATCC", "ACGTTGCAGGTACCTAGG"]
    for order in (reads, reads[::-1]):
        counts = dict(Counter(brute_kmers(order, 3, False)))  # insertion order follows the reads
        got = extract_contigs(build_assembly_graph(counts), 3, both_strands=False)
        assert [(c.id, c.length, c.sequence) for c in got] == expected


_CYCLES_CHILD = """
import json
from collections import Counter
from genomes import brute_kmers
from eulerase.assembly import build_assembly_graph, extract_contigs

counts = dict(Counter(brute_kmers(["ACGTACGT", "AACCAACC"], 4, False)))
got = extract_contigs(build_assembly_graph(counts), 4, both_strands=False)
print(json.dumps([(c.id, c.length, c.sequence) for c in got]))
"""


def test_equal_cycles_are_numbered_the_same_in_two_processes():
    # Each read is a 4-edge cycle written once as its circle: ACGTACGT is ACGT and AACCAACC
    # is AACC; both have length 4, so the tie goes by sequence: AACC is contig 1. Seeds 1 and
    # 4 gave opposite orders before the fix.
    expected = [(1, 4, "AACC"), (2, 4, "ACGT")]
    for seed in ("1", "4"):
        env = child_env() | {
            "PYTHONHASHSEED": seed,
            "PYTHONPATH": os.pathsep.join([child_env()["PYTHONPATH"], str(Path(__file__).parent)]),
        }
        done = subprocess.run(
            [sys.executable, "-c", _CYCLES_CHILD],
            env=env,
            capture_output=True,
            text=True,
            timeout=60,
            check=True,
        )
        assert [tuple(row) for row in json.loads(done.stdout)] == expected


def test_unpaired_twin_is_reported_as_it_is():
    forward = count_kmers(brute_kmers(["ACGTT"], 3, False))
    alone = extract_contigs(build_assembly_graph(forward), 3, both_strands=True)
    assert [(c.sequence, c.strand) for c in alone] == [("ACGTT", "unpaired")]


def test_a_contig_and_its_twin_are_kept_once_and_a_palindrome_is_kept_once():
    paired = count_kmers(brute_kmers(["ACCGT"], 4, True))
    both = extract_contigs(build_assembly_graph(paired), 4, both_strands=True)
    assert [(c.sequence, c.strand) for c in both] == [("ACCGT", "smaller-of-pair")]
    own = count_kmers(brute_kmers(["ATGCAT"], 6, True))
    pal = extract_contigs(build_assembly_graph(own), 6, both_strands=True)
    assert [(c.sequence, c.strand) for c in pal] == [("ATGCAT", "smaller-of-pair")]
    # A cycle that is its own reverse complement: the circle ACGT is ACGT again when reversed
    # and complemented. Its k-mers are ACG, CGT, GTA, TAC, and the reverse strand gives the
    # same four.
    ring = count_kmers(brute_kmers(["ACGTAC"], 3, True))
    cyc = extract_contigs(build_assembly_graph(ring), 3, both_strands=True)
    assert [(c.sequence, c.circular, c.strand) for c in cyc] == [("ACGT", True, "smaller-of-pair")]
    # A cycle that is its own reverse complement only up to rotation. GATCGA at k=3 has
    # k-mers GAT, ATC, TCG, CGA; their reverse complements ATC, GAT, CGA, TCG are the same
    # four. Nodes GA, AT, TC, CG form one cycle whose circle is ATCG (its least rotation);
    # the reverse complement of ATCG is CGAT, a rotation of ATCG but not equal to it.
    # TCGATCG at k=4 has TCGA, CGAT, GATC, ATCG; reverse complements TCGA, ATCG, GATC, CGAT
    # are the same four, and the circle is ATCG again (rc CGAT).
    for read, k, expected in (("GATCGA", 3, "ATCG"), ("TCGATCG", 4, "ATCG")):
        twisted = count_kmers(brute_kmers([read], k, True))
        got = extract_contigs(build_assembly_graph(twisted), k, both_strands=True)
        assert [(c.sequence, c.circular, c.strand) for c in got] == [
            (expected, True, "smaller-of-pair")
        ]


# --- a cycle is its circle, written once ------------------------------------------------------


@pytest.mark.parametrize("both", [False, True], ids=["one-strand", "both-strands"])
@pytest.mark.parametrize(
    ("read", "k", "circle"),
    [
        ("AAAAAAAA", 3, "A"),  # a homopolymer: one letter
        ("ATGCATGCATGC", 3, "ATGC"),  # a tandem repeat: shorter than the input
        ("ACACACACAC", 5, "AC"),  # n = 2, below k - 1 = 4
        ("ACGTAC", 3, "ACGT"),  # its own reverse complement
        ("GATCGA", 3, "ATCG"),  # its reverse complement CGAT is another rotation
        ("TCGATCG", 4, "ATCG"),
    ],
)
def test_a_cycle_is_its_circle_written_once(read, k, circle, both):
    # the letter count must not be mistakable for the read's, nor for the k-1 letters a cycle
    # used to repeat at its end
    assert len(circle) not in (len(read), k - 1)
    (contig,) = contigs_of([read], k, both)
    assert (contig.sequence, contig.circular) == (circle, True)
    assert contig.length == contig.edges == len(circle)


def _expected_circle(unit: str, both: bool) -> str:
    """The least rotation of the circle, or with both strands the smaller of that and the least
    rotation of its reverse complement."""
    least = least_rotation(unit)
    return min(least, least_rotation(revcomp(unit))) if both else least


def _distinct_circle(n: int, k: int, both: bool, seed: int) -> str:
    """n letters whose n cyclic (k-1)-mers are all different (and, with both strands, differ from
    their reverse complements), so that the reads are one cycle of n k-mers."""
    for attempt in range(500):
        unit = random_dna(n, seed * 1000 + attempt)
        ring = unit * (k // n + 2)
        nodes = {ring[i : i + k - 1] for i in range(n)}
        if len(nodes) == n and not (both and nodes & {revcomp(x) for x in nodes}):
            return unit
    raise AssertionError(f"no circle of {n} letters for k={k}, both strands {both}")


@pytest.mark.parametrize("seed", range(120))
def test_a_cycle_matches_a_brute_force_circle(seed):
    rng = random.Random(seed)
    n, both = rng.randint(1, 20), bool(rng.getrandbits(1))
    # k - 1 letters give 4 ** (k - 1) different (k-1)-mers: room for n of them, with spare
    k = rng.randint(next(k for k in range(2, 9) if 4 ** (k - 1) >= 4 * n), 8)
    unit = _distinct_circle(n, k, both, seed)
    read = unit * ((n + k - 1) // n + 1)  # the circle written round, long enough for every k-mer
    contigs = contigs_of([read], k, both)
    assert [(c.sequence, c.circular, c.length, c.edges) for c in contigs] == [
        (_expected_circle(unit, both), True, n, n)
    ]


def test_oracle_rejects_a_rotation_that_is_not_the_least():
    assert least_rotation("TGCA") == "ATGC" != "TGCA"
    assert _expected_circle("ATGC", False) != "CATG"
    # TTTA: least rotation ATTT; its reverse complement TAAA has the smaller least rotation AAAT
    assert _expected_circle("TTTA", False) == "ATTT"
    assert _expected_circle("TTTA", True) == "AAAT"


# --- the metrics count a cycle by its circle, and the mean is over k-mers ----------------------


def _direct_mean(contigs, counts: dict[str, int], k: int) -> float:
    """Times seen over distinct k-mers, from the k-mers the contigs spell and the counts."""
    seen = [counts[w] for c in contigs for w in _windows(c.sequence, k, c.circular)]
    return sum(seen) / len(seen)


def test_metrics_count_a_cycle_by_its_circle():
    reads = ["ATGCATGCATGC", "AACCGGTT"]  # no 2-mer shared: a cycle ATGC and a line AACCGGTT
    contigs = contigs_of(reads, 3, False)
    m = assembly_metrics(contigs)
    assert [(c.sequence, c.circular, c.edges, c.kmer_count) for c in contigs] == [
        ("AACCGGTT", False, 6, 6),
        ("ATGC", True, 4, 10),
    ]
    # the quantities that could be swapped are all different
    assert len({m.total_length, m.longest, m.shortest, m.mean_length}) == 4
    assert (m.total_length, m.longest, m.shortest, m.mean_length) == (12, 8, 4, 6.0)
    assert (m.n50, m.l50, m.n90) == (8, 1, 4)
    assert m.gc == 6 / 12  # G and C over the 12 letters of both contigs
    # 16 k-mers seen over 10 distinct k-mers, not the mean of the contig means
    assert m.mean_multiplicity == 16 / 10
    assert m.mean_multiplicity != statistics.mean(c.mean_multiplicity for c in contigs) == 1.75


def _mean_case(seed: int):
    rng = random.Random(seed)
    alphabet = ("AC", "ACG", "ACGT")[seed % 3]
    k = rng.randint(2, 5)
    reads = []
    for _ in range(rng.randint(1, 4)):
        chunk = "".join(rng.choice(alphabet) for _ in range(rng.randint(k, k + 10)))
        reads.append(chunk * rng.choice((1, 3)))  # three turns of a chunk can close a cycle
    return reads, k, bool(rng.getrandbits(1)), rng.randint(1, 3)


def test_mean_multiplicity_is_total_kmer_count_over_total_edges():
    cycles = differing = 0
    for seed in range(150):
        reads, k, both, minimum = _mean_case(seed)
        counts = count_kmers(brute_kmers(reads, k, both))
        ag, _ = filter_min_multiplicity(build_assembly_graph(counts), minimum)
        contigs = extract_contigs(ag, k, both_strands=both)
        if not contigs:
            continue
        mean = assembly_metrics(contigs).mean_multiplicity
        assert mean == pytest.approx(_direct_mean(contigs, counts, k))
        cycles += any(c.circular for c in contigs)
        differing += mean != pytest.approx(statistics.mean(c.mean_multiplicity for c in contigs))
    # not vacuous: the inputs hold cycles, and the two ways of averaging differ on some
    assert cycles >= 10 and differing >= 10


# --- the links of a cycle and of a contig whose ends overlap ------------------------------------


@pytest.mark.parametrize("both", [False, True], ids=["one-strand", "both-strands"])
@pytest.mark.parametrize(
    ("read", "k"),
    [
        ("ACACACACAC", 5),  # n = 2 is below k - 1 = 4
        ("ATGCATGCATGC", 3),  # n = 4 is above k - 1 = 2
        ("ACGTAC", 3),  # its own reverse complement
    ],
)
def test_a_cycle_has_one_self_link_and_no_other(read, k, both):
    contigs = contigs_of([read], k, both)
    assert [c.circular for c in contigs] == [True]
    assert gfa_links(contigs, k, both) == [(1, "+", 1, "+")]


def test_a_contig_whose_ends_overlap_links_to_itself():
    contigs = contigs_of(["AACGTAA", "AAG"], 3, False)
    assert [(c.id, c.sequence, c.circular) for c in contigs] == [
        (1, "AACGTAA", False),
        (2, "AAG", False),
    ]
    seqs = {c.id: c.sequence for c in contigs}
    # every pair of contigs whose last 2 letters equal the other's first 2, tried one by one
    pairs = {(a, "+", b, "+") for a, s in seqs.items() for b, t in seqs.items() if s[-2:] == t[:2]}
    assert pairs == {(1, "+", 1, "+"), (1, "+", 2, "+")}
    assert set(gfa_links(contigs, 3, False)) == pairs

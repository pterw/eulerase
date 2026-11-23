"""Tests for the assembly graph, contigs, metrics and the k-mer spectrum."""

import statistics
from collections import Counter

import pytest
from genomes import brute_kmers, random_dna, revcomp

from eulerase.assembly import (
    assembly_metrics,
    build_assembly_graph,
    extract_contigs,
    kmer_spectrum,
)
from eulerase.kmers import count_kmers


def _contigs(reads: list[str], k: int, both: bool):
    counts = count_kmers(brute_kmers(reads, k, both))
    return extract_contigs(build_assembly_graph(counts), k, both_strands=both)


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


def _is_the_genome(sequences: list[str], genome: str, both: bool) -> bool:
    allowed = {genome, revcomp(genome)} if both else {genome}
    return len(sequences) == 1 and sequences[0] in allowed


@pytest.mark.parametrize("both", [False, True])
@pytest.mark.parametrize("seed", range(20))
def test_error_free_genome_gives_one_contig_equal_to_it(seed, both):
    k = 11 + seed % 15
    genome = _repeat_free_genome(300 + (seed * 37) % 701, k, both, seed)
    contigs = _contigs(_tiles(genome), k, both)
    assert _is_the_genome([c.sequence for c in contigs], genome, both)
    assert contigs[0].strand == ("smaller-of-pair" if both else "as-assembled")


def test_oracle_rejects_a_contig_one_base_off():
    genome = random_dna(300, seed=1)
    off = genome[:100] + ("A" if genome[100] != "A" else "C") + genome[101:]
    assert _is_the_genome([genome], genome, False)
    assert not _is_the_genome([off], genome, False)
    assert not _is_the_genome([off], genome, True)


@pytest.mark.parametrize(
    ("genome", "k"),
    [
        (random_dna(3000, seed=11), 21),
        (random_dna(400, seed=5), 15),
        ("CTTCGCCTGATACG", 5),
    ],
)
def test_both_strands_gives_one_canonical_contig_per_strand_pair(genome, k):
    contigs = _contigs([genome], k, True)
    seqs = [c.sequence for c in contigs]
    assert len(genome) - 5 <= sum(len(s) for s in seqs) <= len(genome) + 5
    for s in seqs:
        assert revcomp(s) not in seqs or revcomp(s) == s
        assert s == min(s, revcomp(s))


def test_forward_input_is_never_deduplicated():
    read = "CTTCGCCTGATACG"
    with_twin = _contigs([read, revcomp(read)], 5, False)
    assert sorted(c.sequence for c in with_twin) == sorted([read, revcomp(read)])
    unrelated = [random_dna(40, seed=21), random_dna(40, seed=22)]
    assert sorted(c.sequence for c in _contigs(unrelated, 9, False)) == sorted(unrelated)


def test_n_does_not_defeat_twin_removal():
    read = "ACGTNACGTTGCAGGATCCATTGACN"
    both = _contigs([read], 5, True)
    seqs = [c.sequence for c in both]
    assert all(revcomp(s) not in seqs or revcomp(s) == s for s in seqs)
    assert any("N" in s for s in seqs)
    both_strand_graph = build_assembly_graph(count_kmers(brute_kmers([read], 5, True)))
    undeduplicated = extract_contigs(both_strand_graph, 5, both_strands=False)
    assert len(both) < len(undeduplicated)
    # has_n follows the letters of each contig, in both directions.
    mixed = _contigs([read, random_dna(30, seed=31)], 7, False)
    assert {c.has_n for c in mixed} == {True, False}
    assert all(c.has_n == ("N" in c.sequence) for c in mixed)


def test_circular_genome_with_both_strands_is_one_circular_contig():
    genome = random_dna(300, seed=9)
    k = 15
    contigs = _contigs([genome + genome[: k - 1]], k, True)
    assert len(contigs) == 1
    assert contigs[0].circular is True
    assert contigs[0].edges == len(genome)


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


def test_oracle_rejects_a_wrong_n50():
    assert _n_statistic([10, 10, 10, 1], 0.5) == (10, 2)
    assert _n_statistic([10, 10, 10, 1], 0.5) != (1, 4)


@pytest.mark.parametrize("seed", range(150))
def test_metrics_match_direct_computation(seed):
    k = 4 + seed % 12
    reads = [random_dna(k + 5 + (seed * (i + 3)) % 40, seed * 7 + i) for i in range(1 + seed % 5)]
    counts = count_kmers(brute_kmers(reads, k, False) + brute_kmers(reads[:1], k, False))
    contigs = extract_contigs(build_assembly_graph(counts), k, both_strands=False)
    per_contig = []
    for c in contigs:
        seen = [counts[c.sequence[i : i + k]] for i in range(len(c.sequence) - k + 1)]
        assert (c.edges, c.kmer_count, c.min_multiplicity) == (len(seen), sum(seen), min(seen))
        assert c.mean_multiplicity == pytest.approx(sum(seen) / len(seen))
        assert c.gc == (c.sequence.count("G") + c.sequence.count("C")) / len(c.sequence)
        per_contig.append(sum(seen) / len(seen))
    lengths = [len(c.sequence) for c in contigs]
    joined = "".join(c.sequence for c in contigs)
    m = assembly_metrics(contigs)
    assert (m.contigs, m.total_length) == (len(contigs), sum(lengths))
    assert (m.longest, m.shortest) == (max(lengths), min(lengths))
    assert m.mean_length == pytest.approx(sum(lengths) / len(lengths))
    assert (m.n50, m.l50) == _n_statistic(lengths, 0.5)
    assert m.n90 == _n_statistic(lengths, 0.9)[0]
    assert m.gc == pytest.approx((joined.count("G") + joined.count("C")) / len(joined))
    assert m.mean_multiplicity == pytest.approx(statistics.mean(per_contig))


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


def test_worked_case_ten_copies():
    contigs = _contigs(["ATGC"] * 10, 3, False)
    ag = build_assembly_graph(count_kmers(brute_kmers(["ATGC"] * 10, 3, False)))
    assert sorted(d["multiplicity"] for _, _, d in ag.edges(data=True)) == [10, 10]
    assert [c.sequence for c in contigs] == ["ATGC"]


def test_worked_case_repeat_inside_one_sequence():
    read = "ATGCATGCATGC"
    counts = count_kmers(brute_kmers([read], 3, False))
    contigs = extract_contigs(build_assembly_graph(counts), 3, both_strands=False)
    assert Counter(counts.values()) == {3: 2, 2: 2}
    assert [(c.sequence, c.circular, c.edges) for c in contigs] == [("ATGCAT", True, 4)]
    assert contigs[0].mean_multiplicity == 2.5
    assert contigs[0].length < len(read)


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="R1: contigs of equal length are numbered in the order the reads came in",
)
def test_contigs_are_numbered_in_a_stable_order():
    reads = ["ACGTTGCACCATGGATCC", "ACGTTGCAGGTACCTAGG"]
    for order in (reads, reads[::-1]):
        counts = dict(Counter(brute_kmers(order, 3, False)))
        got = extract_contigs(build_assembly_graph(counts), 3, both_strands=False)
        sequences = sorted({c.sequence for c in got}, key=lambda s: (-len(s), s))
        assert [(c.id, c.sequence) for c in got] == list(enumerate(sequences, 1))


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="A5: a contig whose twin is absent is reported as smaller-of-pair",
)
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

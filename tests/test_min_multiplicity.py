"""No k-mer is removed unless the user sets a minimum multiplicity, and the count is honest."""

import random
from collections import Counter

import pytest
from genomes import brute_kmers, is_the_genome, random_dna, revcomp
from pipeline import contig_sequences

from eulerase.analysis import Settings, analyze
from eulerase.sequence import Record

GENOME_LENGTH = 1000


def _genome_without_repeats(k: int, seed: int) -> str:
    """A random genome in which no (k-1)-mer occurs twice, on either strand."""
    for attempt in range(100):
        genome = random_dna(GENOME_LENGTH, seed * 1000 + attempt)
        nodes = brute_kmers([genome], k - 1, True)
        if len(nodes) == len(set(nodes)):
            return genome
    pytest.fail("no genome without a repeated (k-1)-mer found")


def _tiled_reads(genome: str) -> list[str]:
    return [genome[i : i + 100] for i in range(0, GENOME_LENGTH - 100 + 1, 5)]


def _random_reads(genome: str, k: int, seed: int) -> list[str]:
    """30x reads of 60 letters, the first and last on the genome ends, covering every k-mer."""
    wanted = set(brute_kmers([genome], k, False))
    count = 30 * GENOME_LENGTH // 60
    for attempt in range(100):
        rng = random.Random(seed * 1000 + attempt)
        starts = [0, GENOME_LENGTH - 60]
        starts += [rng.randrange(GENOME_LENGTH - 60 + 1) for _ in range(count - 2)]
        reads = [genome[s : s + 60] for s in starts]
        if wanted <= set(brute_kmers(reads, k, False)):
            return reads
    pytest.fail("no read set covering every k-mer found")


@pytest.mark.parametrize("both_strands", [False, True])
@pytest.mark.parametrize("seed", [1, 2, 3])
@pytest.mark.parametrize("k", [15, 21])
@pytest.mark.parametrize("scheme", ["tiled", "random"])
def test_error_free_reads_give_one_contig_equal_to_the_genome(scheme, k, seed, both_strands):
    genome = _genome_without_repeats(k, seed)
    reads = _tiled_reads(genome) if scheme == "tiled" else _random_reads(genome, k, seed)
    contigs = contig_sequences(reads, k, both_strands)
    assert is_the_genome(contigs, genome, both_strands), [len(c) for c in contigs]


@pytest.mark.parametrize("minimum", [2, 3, 4, 5])
def test_removed_count_equals_an_independent_count(minimum):
    k = 15
    rng = random.Random(7)
    genome = random_dna(GENOME_LENGTH, 7)
    reads = []
    for read in _tiled_reads(genome):
        letters = list(read)
        for i in range(len(letters)):
            if rng.random() < 0.01:
                letters[i] = rng.choice([b for b in "ACGT" if b != letters[i]])
        reads.append("".join(letters))
    records = [Record(f"read{i}", read) for i, read in enumerate(reads, 1)]
    seen = Counter(brute_kmers(reads, k, False))

    default = analyze(records, Settings(k=k))
    assert default.removed_kmers == 0
    assert (default.clipped_branches, default.clipped_kmers) == (0, 0)
    assert default.distinct_kmers == len(seen)
    assert analyze(records, Settings(k=k, clip_tips=True)).clipped_branches > 0

    expected = sum(1 for n in seen.values() if n < minimum)
    assert expected > 0
    result = analyze(records, Settings(k=k, min_multiplicity=minimum))
    assert result.removed_kmers == expected


def test_oracle_rejects_a_contig_that_is_short_off_by_a_base_or_the_wrong_strand():
    genome = random_dna(GENOME_LENGTH, 11)
    off = genome[:100] + ("A" if genome[100] != "A" else "C") + genome[101:]
    assert is_the_genome([genome], genome, False)
    assert not is_the_genome([off], genome, False)
    assert not is_the_genome([genome[:-75]], genome, False)
    assert not is_the_genome([revcomp(genome)], genome, False)
    assert not is_the_genome([revcomp(genome)[:-75]], genome, True)

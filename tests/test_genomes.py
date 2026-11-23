"""The random sequence helper that later tests compare against brute-force answers."""

from genomes import random_dna


def test_random_dna_is_a_reproducible_mixed_sequence():
    s = random_dna(200, 1)
    assert len(s) == 200
    assert set(s) == set("ACGT")
    assert random_dna(200, 1) == s
    assert random_dna(200, 2) != s

"""The presenter's pure functions: the k-mer estimate and the window size."""

import random

import pytest
from genomes import brute_kmers, random_dna

from eulerase.gui.presenter import (
    DEFAULT_WINDOW,
    MIN_WINDOW,
    SCREEN_ALLOWANCE,
    estimate_kmers,
    fit_window,
)
from eulerase.sequence import Record


def random_records(seed: int) -> list[Record]:
    rng = random.Random(seed)
    return [Record(f"r{i}", random_dna(rng.randint(1, 40), seed * 10 + i)) for i in range(4)]


@pytest.mark.parametrize("both_strands", [False, True])
@pytest.mark.parametrize("seed", range(6))
def test_the_estimate_equals_a_brute_force_count(seed, both_strands):
    records = random_records(seed)
    k = 2 + seed * 2  # 2 to 12: the longest records are shorter than the largest k
    expected = len(brute_kmers([r.sequence for r in records], k, both_strands))
    assert estimate_kmers(records, k, both_strands=both_strands, input_mode="sequences") == expected


@pytest.mark.parametrize("both_strands", [False, True])
def test_the_estimate_for_a_kmer_list_is_one_per_record(both_strands):
    records = [Record(f"kmer{i}", kmer) for i, kmer in enumerate(["ACG", "CGT", "GTA"], 1)]
    expected = len(records) * (2 if both_strands else 1)
    assert estimate_kmers(records, None, both_strands=both_strands, input_mode="kmers") == expected


def fits(size, screen):
    return all(
        used + taken <= available
        for used, taken, available in zip(size, SCREEN_ALLOWANCE, screen, strict=True)
    )


@pytest.mark.parametrize(
    "screen", [(1366, 768), (1280, 720), (1920, 1080), (800, 600), (1024, 700)]
)
def test_the_window_fits_the_screen_and_never_goes_below_its_minimum(screen):
    size = fit_window(*screen)
    assert all(m <= s <= d for m, s, d in zip(MIN_WINDOW, size, DEFAULT_WINDOW, strict=True))
    if fits(MIN_WINDOW, screen):
        assert fits(size, screen)
    if fits(DEFAULT_WINDOW, screen):
        assert size == DEFAULT_WINDOW
    if not fits(MIN_WINDOW, screen):
        assert size == MIN_WINDOW

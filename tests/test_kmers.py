"""Tests for turning records into k-mers: the k check, the limit, both strands."""

import tracemalloc
from collections import Counter

import pytest
from genomes import brute_kmers, random_dna

from eulerase import kmers
from eulerase.errors import EuleraseError
from eulerase.sequence import Record, parse_records


def _records(*sequences: str) -> list[Record]:
    return [Record(f"r{i}", s) for i, s in enumerate(sequences, 1)]


def _same_multiset(actual: list[str] | tuple[str, ...], expected: list[str]) -> bool:
    return Counter(actual) == Counter(expected)


@pytest.mark.parametrize("unique", [False, True])
@pytest.mark.parametrize("both", [False, True])
@pytest.mark.parametrize("k", range(2, 13))
def test_from_sequences_matches_brute_force_counts(k, both, unique):
    seqs = [random_dna(40 + 7 * k, seed=k), random_dna(25, seed=100 + k), "ACG"]
    result = kmers.from_sequences(_records(*seqs), k, both_strands=both, unique=unique)
    expected = brute_kmers(seqs, k, both)
    if unique:
        assert list(result.items) == sorted(set(expected))
    else:
        assert list(result.items) == expected
    assert result.k == k


def test_oracle_rejects_a_missing_kmer():
    expected = brute_kmers(["ACGTACGT"], 3, False)
    assert _same_multiset(expected, expected)
    assert not _same_multiset(expected[1:], expected)


def test_self_reverse_complement_kmer_counts_twice():
    result = kmers.from_sequences(_records("ATGCAT"), 6, both_strands=True, unique=False)
    assert Counter(result.items)["ATGCAT"] == 2


def test_short_records_are_counted_not_dropped_silently():
    result = kmers.from_sequences(
        _records("ACGTACGT", "AC", "GGG"), 4, both_strands=False, unique=False
    )
    assert result.short_records == 2
    assert list(result.items) == brute_kmers(["ACGTACGT"], 4, False)
    with pytest.raises(EuleraseError, match="k=4"):
        kmers.from_sequences(_records("AC", "GG"), 4, both_strands=False, unique=False)


def test_every_record_empty_is_refused_with_its_own_sentence():
    # No k helps an empty record, so the refusal must not tell the user to use a smaller k.
    with pytest.raises(EuleraseError) as empty:
        kmers.from_sequences(_records("", ""), 3, both_strands=False, unique=False)
    assert "empty" in str(empty.value)
    assert "smaller k" not in str(empty.value)
    # a header with nothing under it is an empty record too, not "no sequence"
    with pytest.raises(EuleraseError, match="Every record is empty"):
        kmers.from_sequences(parse_records(">only a header\n"), 3, both_strands=False, unique=False)
    beside_long = kmers.from_sequences(_records("", "ACGTAC"), 3, both_strands=False, unique=False)
    assert beside_long.short_records == 1
    assert list(beside_long.items) == brute_kmers(["ACGTAC"], 3, False)
    # an empty record beside a record that is only too short: a smaller k does help here
    with pytest.raises(EuleraseError, match="smaller k"):
        kmers.from_sequences(_records("", "AC"), 3, both_strands=False, unique=False)


def test_no_records_is_refused():
    with pytest.raises(EuleraseError):
        kmers.from_sequences([], 3, both_strands=False, unique=False)
    with pytest.raises(EuleraseError):
        kmers.from_kmer_list([], both_strands=False, unique=False)


@pytest.mark.parametrize("bad_k", [1, 0, -2, 2.5, True])
def test_invalid_k_is_refused_everywhere(bad_k):
    with pytest.raises(EuleraseError, match="k"):
        kmers.validate_k(bad_k)
    with pytest.raises(EuleraseError, match="k"):
        kmers.from_sequences(_records("ACGTACGT"), bad_k, both_strands=False, unique=False)


@pytest.mark.parametrize("sequences", [["A", "C"], [""]])
def test_kmer_list_with_an_invalid_length_is_refused(sequences):
    with pytest.raises(EuleraseError, match="k"):
        kmers.from_kmer_list(_records(*sequences), both_strands=False, unique=False)


@pytest.mark.parametrize(
    ("mode", "size", "both", "accepted"),
    [
        ("sequences", 14, False, True),  # 10 k-mers of length 5
        ("sequences", 15, False, False),  # 11
        ("sequences", 8, True, True),  # 5 + 5 with k=4
        ("sequences", 9, True, False),  # 6 + 6
        ("list", 10, False, True),
        ("list", 11, False, False),
        ("list", 5, True, True),  # 5 + 5 reverse complements
        ("list", 6, True, False),
    ],
)
def test_limit_cases(monkeypatch, mode, size, both, accepted):
    monkeypatch.setattr(kmers, "MAX_KMERS", 10)

    def build():
        if mode == "sequences":
            k = 4 if both else 5
            return kmers.from_sequences(
                _records(random_dna(size, seed=1)), k, both_strands=both, unique=False
            )
        words = [random_dna(4, seed=i) for i in range(size)]
        return kmers.from_kmer_list(_records(*words), both_strands=both, unique=False)

    if accepted:
        assert len(build().items) <= 10
    else:
        with pytest.raises(EuleraseError, match="10"):
            build()


def test_oversized_input_is_refused_before_the_list_is_built(monkeypatch):
    monkeypatch.setattr(kmers, "MAX_KMERS", 1000)
    record = _records("ACGT" * 5_000)
    k = 10_000
    count = 2 * (len(record[0].sequence) - k + 1)  # both strands
    peak_bound = 1_000_000
    assert count == 20_002
    assert count != 1_000 + 1  # a message made by stopping early cannot name the count
    assert (1_000 + 1) * k > peak_bound  # building strings until the limit cannot stay under
    tracemalloc.start()
    try:
        with pytest.raises(EuleraseError) as refused:
            kmers.from_sequences(record, k, both_strands=True, unique=False)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert "20,002" in str(refused.value)
    assert "1,000" in str(refused.value)
    assert peak < peak_bound


def test_both_strands_refuses_non_dna_naming_the_letters():
    with pytest.raises(EuleraseError, match="X"):
        kmers.from_sequences(_records("ACGXT"), 3, both_strands=True, unique=False)
    with pytest.raises(EuleraseError, match="X"):
        kmers.from_kmer_list(_records("ACX", "CGT"), both_strands=True, unique=False)
    # Without both strands the same letters are fine.
    kmers.from_sequences(_records("ACGXT"), 3, both_strands=False, unique=False)


def test_kmer_list_mixed_lengths_refused():
    with pytest.raises(EuleraseError, match="same length"):
        kmers.from_kmer_list(_records("ACG", "ACGT"), both_strands=False, unique=False)


def test_kmer_list_uses_its_own_length_and_adds_reverse_complements():
    result = kmers.from_kmer_list(_records("AAC", "GGT"), both_strands=True, unique=False)
    assert result.k == 3
    assert list(result.items) == ["AAC", "GGT", "GTT", "ACC"]


def test_unique_reports_duplicates_dropped():
    result = kmers.from_sequences(_records("AAAAA"), 2, both_strands=False, unique=True)
    assert result.items == ("AA",)
    assert result.duplicates_dropped == 3
    kept = kmers.from_sequences(_records("AAAAA"), 2, both_strands=False, unique=False)
    assert kept.duplicates_dropped == 0


def test_count_kmers_counts_and_sorts_keys():
    items = ["GT", "AC", "GT", "AC", "AA", "GT"]
    counts = kmers.count_kmers(items)
    assert list(counts) == sorted(set(items))
    assert counts == {"AA": 1, "AC": 2, "GT": 3}

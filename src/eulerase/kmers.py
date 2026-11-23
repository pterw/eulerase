"""Turn records into k-mers, with the one check of k and the one input limit."""

from collections import Counter
from collections.abc import Iterator, Sequence
from dataclasses import dataclass

from eulerase.errors import EuleraseError
from eulerase.sequence import Record, require_dna, reverse_complement

MAX_KMERS = 500_000


@dataclass(frozen=True)
class Kmers:
    items: tuple[str, ...]
    k: int
    short_records: int
    duplicates_dropped: int


def validate_k(k: object) -> int:
    """Return k when it is an int (not a bool) of at least 2."""
    if isinstance(k, bool) or not isinstance(k, int) or k < 2:
        raise EuleraseError(f"k must be a whole number of at least 2, not {k!r}; give a larger k.")
    return k


def _limit_error(n: int) -> EuleraseError:
    return EuleraseError(
        f"Input produced more than {MAX_KMERS:,} k-mers (counted {n:,} before stopping); "
        f"the limit is {MAX_KMERS:,}. Use a smaller input or a larger k."
    )


def _iter_kmers(sequences: Sequence[str], k: int, both: bool) -> Iterator[str]:
    """Yield each sequence's k-mers, then its reverse complement's when `both`, lazily so
    the caller can stop as soon as the limit is passed."""
    for seq in sequences:
        for s in (seq, reverse_complement(seq)) if both else (seq,):
            for i in range(len(s) - k + 1):
                yield s[i : i + k]


def _finish(items: list[str], k: int, short_records: int, unique: bool) -> Kmers:
    if not unique:
        return Kmers(tuple(items), k, short_records, 0)
    distinct = sorted(set(items))
    return Kmers(tuple(distinct), k, short_records, len(items) - len(distinct))


def from_sequences(records: Sequence[Record], k: int, *, both_strands: bool, unique: bool) -> Kmers:
    """The k-mers of every record that is at least k long (and of its reverse complement
    when `both_strands`); shorter records are counted in `short_records`."""
    k = validate_k(k)
    if not records:
        raise EuleraseError("There is no sequence to read; give a sequence or a file.")
    sequences = [r.sequence for r in records]
    if both_strands:
        require_dna(sequences)
    long_enough = [s for s in sequences if len(s) >= k]
    if not long_enough:
        raise EuleraseError(f"No record is at least k={k} letters long; use a smaller k.")
    items: list[str] = []
    for km in _iter_kmers(long_enough, k, both_strands):
        items.append(km)
        if len(items) > MAX_KMERS:
            raise _limit_error(len(items))
    return _finish(items, k, len(sequences) - len(long_enough), unique)


def from_kmer_list(records: Sequence[Record], *, both_strands: bool, unique: bool) -> Kmers:
    """Each record's sequence is one k-mer; all must have one length, which is k."""
    if not records:
        raise EuleraseError("There are no k-mers to read; give some k-mers.")
    words = [r.sequence for r in records]
    if len({len(w) for w in words}) != 1:
        raise EuleraseError("All k-mers must have the same length; fix or remove the odd ones.")
    k = validate_k(len(words[0]))
    total = len(words) * (2 if both_strands else 1)
    if total > MAX_KMERS:
        raise _limit_error(total)
    if both_strands:
        require_dna(words)
        words = words + [reverse_complement(w) for w in words]
    return _finish(words, k, 0, unique)


def count_kmers(items: Sequence[str]) -> dict[str, int]:
    """Occurrences of each k-mer, keys in ascending string order."""
    return dict(sorted(Counter(items).items()))

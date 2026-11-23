"""Six-frame ORF finder, checked against an independent per-frame Biopython translation."""

import random

import pytest
from Bio.Seq import Seq
from genomes import revcomp

from eulerase.errors import EuleraseError
from eulerase.orfs import MIN_ORF_AA, Orf, find_orfs, orf_rule_text


def _frame_pieces(s: str, f: int) -> list[tuple[int, str]]:
    """Translate frame f of s whole; return (aa offset, aa text) for each stop-terminated piece."""
    frame = s[f:]
    frame = frame[: len(frame) - len(frame) % 3]
    pieces, offset = [], 0
    for piece in str(Seq(frame).translate()).split("*")[:-1]:
        pieces.append((offset, piece))
        offset += len(piece) + 1
    return pieces


def six_frame_orfs(seq: str, min_aa: int) -> list[Orf]:
    """Oracle: split each translated frame at stops, start each piece at its first M."""
    out = []
    for strand, s in ((1, seq), (-1, revcomp(seq))):
        n = len(s)
        for f in range(3):
            for offset, piece in _frame_pieces(s, f):
                m = piece.find("M")
                if m < 0 or len(piece) - m < min_aa:
                    continue
                aa = piece[m:]
                lo = f + 3 * (offset + m)
                hi = lo + 3 * len(aa) + 3
                start, end = (lo, hi) if strand == 1 else (n - hi, n - lo)
                out.append(Orf(strand, strand * (f + 1), start, end, 3 * len(aa) + 3, len(aa), aa))
    return sorted(out, key=lambda o: (-o.strand, abs(o.frame), o.start))


def _random_with_n(rng: random.Random, n: int) -> str:
    return "".join(rng.choices("ACGTN", weights=(10, 10, 10, 10, 1), k=n))


@pytest.mark.parametrize("min_aa", [10, 3])
def test_find_orfs_matches_six_frame_scan(min_aa):
    """Every field of every ORF equals the per-frame translation, several ATGs per stop included."""
    rng = random.Random(7)
    found = 0
    for _ in range(300):
        seq = _random_with_n(rng, rng.randint(30, 300))
        expected = six_frame_orfs(seq, min_aa)
        assert list(find_orfs(seq, min_aa=min_aa)) == expected, seq
        found += len(expected)
    assert found > 50  # the comparison is not vacuous


_FORWARD = "ATG" + "GCT" * 12 + "TAA"  # 42 nt: ATG, 12 x GCT, TAA
# Reverse case: "CC" + revcomp(_FORWARD), 44 nt. Its reverse complement is _FORWARD + "GG", so the
# ORF sits in the first reading frame of the reverse strand, strand -1, frame -1, at 0..42 there,
# which is 44 - 42 = 2 to 44 - 0 = 44 on the forward strand. Length 42 nt, 13 aa, M + 12 A.
_REVERSE = "CC" + revcomp(_FORWARD)


@pytest.mark.parametrize(
    "seq, expected",
    [
        (_FORWARD, Orf(1, 1, 0, 42, 42, 13, "M" + "A" * 12)),
        (_REVERSE, Orf(-1, -1, 2, 44, 42, 13, "M" + "A" * 12)),
    ],
)
def test_oracle_matches_hand_worked_orfs(seq, expected):
    assert six_frame_orfs(seq, 3) == [expected]


def test_one_orf_per_stop_codon():
    seq = "ATG" + "ATG" + "GCT" * 12 + "TAA"
    orfs = find_orfs(seq, min_aa=3)
    assert len(orfs) == 1
    assert orfs[0].start == 0 and orfs[0].length_aa == 14


def test_orf_running_off_the_end_is_not_reported():
    assert find_orfs("ATG" + "GCT" * 20, min_aa=3) == ()
    assert len(find_orfs("ATG" + "GCT" * 20 + "TAA", min_aa=3)) == 1


def test_non_dna_letters_are_refused():
    with pytest.raises(EuleraseError):
        find_orfs("ATGQQQTAA")


def test_orf_rule_text_states_start_codon_code_and_minimum():
    text = orf_rule_text()
    assert "ATG" in text and "standard" in text and str(MIN_ORF_AA) in text

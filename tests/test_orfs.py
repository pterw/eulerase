"""Six-frame ORF finder, checked against an independent per-frame Biopython translation and,
for circles and the N rule, a brute-force reference."""

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


@pytest.mark.parametrize("min_aa", [10, 3])
def test_find_orfs_matches_six_frame_scan(min_aa):
    """Every field of every ORF equals the per-frame translation, several ATGs per stop included."""
    rng = random.Random(7)
    found = 0
    for _ in range(300):
        seq = "".join(rng.choices("ACGT", k=rng.randint(30, 300)))
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


_STOP_CODONS = ("TAA", "TAG", "TGA")


def brute_orfs(sequence: str, circular: bool, min_aa: int) -> list[Orf]:
    """Reference: from every ATG read codons to the first stop, one longest ORF per stop codon."""
    n = len(sequence)
    found = []
    for strand, s in ((1, sequence), (-1, revcomp(sequence))):
        longest = {}  # stop codon position (on the circle, when circular) -> (first letter, length)
        for p in range(n):
            if "".join(s[(p + i) % n] for i in range(3)) != "ATG" or (not circular and p + 3 > n):
                continue
            q = p
            while True:
                reach = q - p + 3  # letters used if the codon at q is read
                if reach > (n if circular else n - p):
                    break  # runs off the end of a line, or would use a letter of a circle twice
                codon = "".join(s[(q + i) % n] for i in range(3))
                if "N" in codon:
                    break
                if codon in _STOP_CODONS:
                    key = q % n if circular else q
                    if key not in longest or reach > longest[key][1]:
                        longest[key] = (p, reach)
                    break
                q += 3
        for p, length in longest.values():
            nt = "".join(s[(p + i) % n] for i in range(length))
            protein = str(Seq(nt[:-3]).translate())
            if len(protein) < min_aa:
                continue
            start = p if strand == 1 else (n - p - length) % n
            found.append(
                Orf(
                    strand,
                    strand * (p % 3 + 1),
                    start,
                    start + length,
                    length,
                    len(protein),
                    protein,
                )
            )
    return sorted(found, key=lambda o: (-o.strand, abs(o.frame), o.start))


_SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"]
_SENSE = [c for c in _SENSE if c not in _STOP_CODONS]


def _circle_sample() -> list[tuple[str, bool]]:
    """Random ACGTN sequences, and planted genes rotated to every cut, as circles and as lines."""
    rng = random.Random(2024)
    sample = []
    for circular in (True, False):
        for _ in range(300):
            letters = rng.choices("ACGTN", weights=(10, 10, 10, 10, 1), k=rng.randint(1, 60))
            sample.append(("".join(letters), circular))
        for _ in range(25):
            gene = (
                "ATG"
                + "".join(rng.choices(_SENSE, k=rng.randint(1, 15)))
                + rng.choice(_STOP_CODONS)
            )
            circle = gene + "".join(rng.choices("ACGT", k=rng.randint(0, 9)))
            if rng.random() < 0.3:  # an N somewhere, to meet the N rule inside a planted gene
                i = rng.randrange(len(circle))
                circle = circle[:i] + "N" + circle[i + 1 :]
            sample.extend((circle[r:] + circle[:r], circular) for r in range(len(circle)))
    return sample


def test_find_orfs_matches_a_brute_force_reference_on_circles_and_lines():
    """Every field of every ORF equals the reference, crossing position 1 and N codons included."""
    sample = _circle_sample()
    crossing = whole = with_n = 0
    seen_strands, seen_frames = set(), set()
    for seq, circular in sample:
        expected = brute_orfs(seq, circular, 1)
        assert list(find_orfs(seq, circular=circular, min_aa=1)) == expected, (seq, circular)
        crossing += sum(o.end > len(seq) for o in expected)
        whole += sum(o.length_nt == len(seq) for o in expected)
        with_n += len(expected) if "N" in seq else 0
        seen_strands |= {o.strand for o in expected}
        seen_frames |= {o.frame for o in expected}
    # the comparison is not vacuous: both strands, all six frames, crossings, full circles, N
    assert seen_strands == {1, -1} and seen_frames == {1, 2, 3, -1, -2, -3}
    assert crossing > 0 and whole > 0 and with_n > 0


_GENE = "ATG" + "GCT" * 6 + "TAA"  # 24 letters, 7 amino acids
_GENE_26 = _GENE + "CC"  # 26 letters: n is not a multiple of 3
_P = "M" + "A" * 6


def _rotate(s: str, r: int) -> str:
    return s[r:] + s[:r]


@pytest.mark.parametrize(
    "seq, expected",
    [
        (_GENE, (Orf(1, 1, 0, 24, 24, 7, _P),)),  # a full circle: end equals the length
        (_rotate(_GENE, 5), (Orf(1, 2, 19, 43, 24, 7, _P),)),  # crosses position 1: end 43 > 24
        (revcomp(_GENE), (Orf(-1, -1, 0, 24, 24, 7, _P),)),
        (_rotate(revcomp(_GENE), 7), (Orf(-1, -2, 17, 41, 24, 7, _P),)),  # crossing, minus strand
        (_GENE_26, (Orf(1, 1, 0, 24, 24, 7, _P),)),
        (_rotate(_GENE_26, 4), (Orf(1, 2, 22, 46, 24, 7, _P),)),
        ("ATGTAA", (Orf(1, 1, 0, 6, 6, 1, "M"),)),  # the smallest circle with an ORF
        ("A", ()),
        ("AT", ()),
        ("ATGTA", ()),
        ("ATG" + "GCT" * 7, ()),  # no stop codon in a whole turn
    ],
)
def test_cycle_cases_worked_by_hand(seq, expected):
    for orf in expected:
        assert orf.start != orf.end and orf.start != abs(orf.frame)  # start, end, frame swappable
    assert find_orfs(seq, circular=True, min_aa=1) == expected


def test_every_rotation_of_a_circle_gives_its_one_orf():
    """A stop codon straddling the cut is still found, once."""
    for r in range(len(_GENE)):
        orfs = find_orfs(_rotate(_GENE, r), circular=True, min_aa=1)
        assert [(o.strand, o.length_nt, o.protein) for o in orfs] == [(1, 24, _P)], r


@pytest.mark.parametrize(
    "seq, expected",
    [
        ("ATGNAAATG" + "GCT" * 3 + "TAA", (Orf(1, 1, 6, 21, 15, 4, "MAAA"),)),  # N ends the first
        ("NNNNNNNNN", ()),
    ],
)
def test_a_codon_with_n_ends_an_orf_on_a_line(seq, expected):
    assert find_orfs(seq, min_aa=1) == expected


def _plain_double_read(seq: str, circular: bool, min_aa: int) -> list[Orf]:
    """Wrong on purpose: the linear finder on s + s, keeping ORFs that start in the first copy."""
    found = find_orfs(seq + seq, min_aa=min_aa) if circular else find_orfs(seq, min_aa=min_aa)
    return [o for o in found if o.start < len(seq)]


def _read_through_n(seq: str, circular: bool, min_aa: int) -> list[Orf]:
    """Wrong on purpose: an N codon is read as an amino acid (N as C is never a stop or an ATG)."""
    return list(find_orfs(seq.replace("N", "C"), circular=circular, min_aa=min_aa))


@pytest.mark.parametrize("wrong", [_plain_double_read, _read_through_n])
def test_oracle_rejects_a_plain_double_read_and_a_read_through_n(wrong):
    sample = _circle_sample()
    differing = sum(wrong(seq, c, 1) != brute_orfs(seq, c, 1) for seq, c in sample)
    assert differing > 0, f"the reference accepted {wrong.__name__} on all {len(sample)} sequences"

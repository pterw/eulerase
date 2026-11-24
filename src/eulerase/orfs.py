"""Six-frame open reading frames: ATG to the first in-frame stop, standard genetic code."""

from collections import deque
from dataclasses import dataclass

from Bio.Seq import Seq

from eulerase.errors import EuleraseError

MIN_ORF_AA = 25
_STOPS = ("TAA", "TAG", "TGA")


@dataclass(frozen=True)
class Orf:
    strand: int
    frame: int
    start: int
    end: int
    length_nt: int
    length_aa: int
    protein: str


def orf_rule_text() -> str:
    return (
        f"ORFs: start codon ATG, standard genetic code, at least {MIN_ORF_AA} amino acids, "
        "the longest per stop codon, ending in a stop codon (none read off an end or across a "
        "codon with N); an ORF is a reading frame, not evidence of a gene."
    )


def _scan_strand(s: str, strand: int, circular: bool, min_aa: int) -> list[Orf]:
    """ORFs of one strand. A circle is read over 2n + 2 letters and only the stop codons of the
    second turn are reported, each of which then sees a whole turn behind it."""
    n = len(s)
    text = (s * 3)[: 2 * n + 2] if circular else s
    first_stop = n if circular else 0
    orfs = []
    for frame in range(3):
        atgs = deque()  # ATG positions since the last stop or N codon in this frame
        for p in range(frame, len(text) - 2, 3):
            codon = text[p : p + 3]
            if codon == "ATG":
                atgs.append(p)
            elif "N" in codon:
                atgs.clear()
            elif codon in _STOPS:
                while atgs and p + 3 - atgs[0] > n:  # an ORF never uses a letter twice
                    atgs.popleft()
                if atgs and p >= first_stop and (p - atgs[0]) // 3 >= min_aa:
                    orfs.append(_make_orf(text[atgs[0] : p + 3], atgs[0] % n, strand, n))
                atgs.clear()
    return orfs


def _make_orf(nt: str, first: int, strand: int, n: int) -> Orf:
    """The ORF whose first letter is at `first` on its own strand; `end` is not reduced modulo n."""
    protein = str(Seq(nt).translate(to_stop=True))
    start = first if strand == +1 else (n - first - len(nt)) % n
    return Orf(
        strand, strand * (first % 3 + 1), start, start + len(nt), len(nt), len(protein), protein
    )


def find_orfs(
    sequence: str, *, circular: bool = False, min_aa: int = MIN_ORF_AA
) -> tuple[Orf, ...]:
    seq_u = sequence.upper()
    if not set(seq_u).issubset("ACGTN"):
        raise EuleraseError(
            "The sequence has letters other than A, C, G, T and N; give DNA to search for ORFs."
        )
    strands = [(seq_u, +1), (str(Seq(seq_u).reverse_complement()), -1)]
    orfs = [o for s, strand in strands for o in _scan_strand(s, strand, circular, min_aa)]
    return tuple(sorted(orfs, key=lambda o: (-o.strand, abs(o.frame), o.start)))

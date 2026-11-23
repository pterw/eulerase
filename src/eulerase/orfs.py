"""Six-frame open reading frames: ATG to the first in-frame stop, standard genetic code."""

from dataclasses import dataclass

from Bio.Seq import Seq

from eulerase.errors import EuleraseError

MIN_ORF_AA = 10
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
        f"ORFs start at ATG, use the standard genetic code, end at a stop codon "
        f"and are at least {MIN_ORF_AA} amino acids long."
    )


def find_orfs(sequence: str, *, min_aa: int = MIN_ORF_AA) -> tuple[Orf, ...]:
    seq_u = sequence.upper()
    if not set(seq_u).issubset("ACGTN"):
        raise EuleraseError(
            "The sequence has letters other than A, C, G, T and N; give DNA to search for ORFs."
        )
    orfs = []
    strands = [(seq_u, +1), (str(Seq(seq_u).reverse_complement()), -1)]
    for s_seq, strand in strands:
        n = len(s_seq)
        for frame in range(3):
            resume = 0  # first position after the stop that ended the previous ORF in this frame
            for i in range(frame, n - 2, 3):
                if i < resume or s_seq[i : i + 3] != "ATG":
                    continue
                for j in range(i + 3, n - 2, 3):
                    if s_seq[j : j + 3] in _STOPS:
                        resume = j + 3
                        aa_len = (j - i) // 3
                        if aa_len >= min_aa:
                            nt_sub = s_seq[i : j + 3]
                            aa_seq = str(Seq(nt_sub).translate(to_stop=True))
                            if strand == +1:
                                start_pos, end_pos = i, j + 3
                            else:
                                start_pos, end_pos = n - (j + 3), n - i
                            orfs.append(
                                Orf(
                                    strand,
                                    strand * (frame + 1),
                                    start_pos,
                                    end_pos,
                                    len(nt_sub),
                                    len(aa_seq),
                                    aa_seq,
                                )
                            )
                        break
                else:
                    break  # no stop before the end: no later ATG in this frame has one either
    return tuple(sorted(orfs, key=lambda o: (-o.strand, abs(o.frame), o.start)))

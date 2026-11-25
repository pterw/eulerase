"""Sequence helpers: records, parsing, alphabet check, reverse complement, GC content."""

import re
from dataclasses import dataclass

from Bio.Seq import Seq

from eulerase.errors import EuleraseError


@dataclass(frozen=True)
class Record:
    name: str
    sequence: str


def reverse_complement(sequence: str) -> str:
    """Reverse complement via Biopython; N stays N."""
    return str(Seq(sequence).reverse_complement())


def parse_records(text: str) -> list[Record]:
    """FASTA-aware parser: with '>' headers, wrapped lines join into one record;
    without headers, every non-empty line is its own record named seq1, seq2, ..."""
    lines = [ln.strip() for ln in text.splitlines()]
    if any(ln.startswith(">") for ln in lines):
        blocks: list[tuple[str, list[str]]] = []
        for ln in lines:
            if not ln:
                continue
            if ln.startswith(">"):
                blocks.append((ln[1:].strip(), []))
            else:
                if not blocks:
                    blocks.append(("", []))
                blocks[-1][1].append(ln)
        found = [(name, "".join(parts)) for name, parts in blocks]
    else:
        found = [("", ln) for ln in lines if ln]
    return [
        Record(name or f"seq{i}", re.sub(r"\s+", "", seq).upper())
        for i, (name, seq) in enumerate(found, 1)
    ]


def parse_kmer_list(text: str) -> list[Record]:
    """Drop '>' lines whole, split the rest on whitespace, commas and semicolons."""
    body = "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith(">"))
    tokens = [t.upper() for t in re.split(r"[\s,;]+", body.strip()) if t]
    return [Record(f"kmer{i}", t) for i, t in enumerate(tokens, 1)]


def is_dna(sequence: str) -> bool:
    """True for a non-empty sequence whose letters are all in ACGTN."""
    return bool(sequence) and set(sequence.upper()) <= set("ACGTN")


def gc_content(sequence: str) -> float | None:
    """(G + C) / (A + C + G + T), N left out; None unless DNA with a known base."""
    if not is_dna(sequence):
        return None
    upper = sequence.upper()
    known = sum(upper.count(base) for base in "ACGT")
    return (upper.count("G") + upper.count("C")) / known if known else None


def require_dna(items: list[str]) -> None:
    """Raise for any item with a letter outside ACGTN (both strands need DNA)."""
    for x in items:
        bad = set(x.upper()) - set("ACGTN")
        if bad:
            raise EuleraseError(
                f"Both strands only makes sense for DNA, but found {sorted(bad)} in '{x}'; "
                "turn both strands off or give DNA."
            )

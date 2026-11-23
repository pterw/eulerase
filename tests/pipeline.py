"""Test support: the one place the tests run reads through the analysis pipeline."""

from eulerase.analysis import Settings, analyze
from eulerase.sequence import Record


def contig_sequences(reads: list[str], k: int, both_strands: bool) -> list[str]:
    records = [Record(f"read{i}", read) for i, read in enumerate(reads, 1)]
    analysis = analyze(records, Settings(k=k, both_strands=both_strands))
    return [contig.sequence for contig in analysis.contigs]

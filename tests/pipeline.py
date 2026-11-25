"""Test support: the one place the tests run reads through the analysis pipeline."""

from genomes import brute_kmers

from eulerase.analysis import Settings, analyze
from eulerase.assembly import Contig, build_assembly_graph, extract_contigs
from eulerase.kmers import count_kmers
from eulerase.sequence import Record


def contigs_of(reads: list[str], k: int, both: bool) -> tuple[Contig, ...]:
    """The contigs of the reads' k-mers, built straight from the assembly graph."""
    counts = count_kmers(brute_kmers(reads, k, both))
    return extract_contigs(build_assembly_graph(counts), k, both_strands=both)


def contig_sequences(reads: list[str], k: int, both_strands: bool) -> list[str]:
    records = [Record(f"read{i}", read) for i, read in enumerate(reads, 1)]
    analysis = analyze(records, Settings(k=k, both_strands=both_strands))
    return [contig.sequence for contig in analysis.contigs]

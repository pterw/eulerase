"""Reading and writing files: sequence file reader (FASTA, FASTQ, GenBank, EMBL; plain or gzip)."""

import gzip
import io
import os
import re
import zlib

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from eulerase.analysis import Analysis, format_cell, settings_text
from eulerase.errors import EuleraseError
from eulerase.orfs import orf_rule_text
from eulerase.sequence import Record, parse_records

_GZIP_MAGIC = b"\x1f\x8b"
_BOM = "\ufeff"
_FIRST_LINE_FORMATS = (("LOCUS", "genbank"), ("ID ", "embl"), (">", "fasta"), ("@", "fastq"))


def _text_of(data: bytes) -> str:
    """Decode, then drop the byte-order mark and any leading blank lines (once, for all users)."""
    text = data.decode("utf-8", errors="replace").removeprefix(_BOM)
    return re.sub(r"\A(?:[^\S\n]*\n)+", "", text)


def detect_format(data: bytes) -> str:
    """Name the format of (already decompressed) file content from its first non-blank line."""
    text = _text_of(data)
    first = next((ln for ln in text.splitlines() if ln.strip()), None)
    if first is None:
        raise EuleraseError("The file is empty; give a file that holds sequences.")
    for prefix, fmt in _FIRST_LINE_FORMATS:
        if first.startswith(prefix):
            return fmt
    if re.fullmatch(r"[A-Za-z*\-\s]+", text):
        return "bare"
    raise EuleraseError(
        "The file is not FASTA, FASTQ, GenBank, EMBL or a bare sequence; check that it is a "
        "sequence file."
    )


def read_records(path: str | os.PathLike) -> list[Record]:
    """Read FASTA / FASTQ / GenBank / EMBL / a bare sequence, plain or gzip (by content)."""
    name = os.fspath(path)
    if not os.path.exists(name):
        raise EuleraseError(f"The file '{name}' does not exist; check the path.")
    try:
        with open(name, "rb") as handle:
            data = handle.read()
    except OSError as err:
        raise EuleraseError(
            f"The file '{name}' could not be read ({err.strerror or 'unreadable'}); "
            "check the path and permissions."
        ) from None
    if data.startswith(_GZIP_MAGIC):
        try:
            data = gzip.decompress(data)
        except (OSError, EOFError, zlib.error):
            raise EuleraseError(
                f"The file '{name}' is damaged or incomplete; download or copy it again."
            ) from None
    fmt = detect_format(data)
    text = _text_of(data)
    if fmt == "bare":
        found = parse_records(text)
    else:
        try:
            parsed = SeqIO.parse(io.StringIO(text), fmt)
            found = [Record(rec.id, str(rec.seq).upper()) for rec in parsed]
        except ValueError:
            raise EuleraseError(
                f"The file '{name}' could not be read as {fmt.upper()}; "
                "check that it is a valid file."
            ) from None
    if not found:
        raise EuleraseError(f"No usable {fmt.upper()} records found in '{name}'.")
    return found


def _seqrecords(a: Analysis) -> list[SeqRecord]:
    """One Biopython record per contig; the settings and, with ORFs, the ORF rule travel in the
    description (FASTA, GenBank) and in the GenBank comment."""
    settings = settings_text(a.settings)
    rule = [] if a.orfs is None else [orf_rule_text()]
    found = [None] * len(a.contigs) if a.orfs is None else a.orfs
    records = []
    for contig, orf_list in zip(a.contigs, found, strict=True):
        count = format_cell(None if orf_list is None else len(orf_list))
        description = " ".join(
            [
                f"len={format_cell(contig.length)} cov={format_cell(contig.mean_multiplicity)}x",
                f"gc={format_cell(contig.gc)} orfs={count}",
                settings,
                *rule,
            ]
        )
        name = f"ctg{contig.id}"
        rec = SeqRecord(Seq(contig.sequence), id=name, name=name, description=description)
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["date"] = "01-JAN-1980"
        rec.annotations["topology"] = "circular" if contig.circular else "linear"
        rec.annotations["comment"] = "\n".join([settings, *rule])
        for orf in orf_list or ():
            rec.features.append(
                SeqFeature(
                    SimpleLocation(orf.start, orf.end, strand=orf.strand),
                    type="CDS",
                    qualifiers={
                        "translation": [orf.protein],
                        "product": ["predicted open reading frame"],
                        "frame": [str(orf.frame)],
                        "codon_start": [1],
                    },
                )
            )
        records.append(rec)
    return records


def _written(a: Analysis, fmt: str) -> str:
    out = io.StringIO()
    SeqIO.write(_seqrecords(a), out, fmt)
    return out.getvalue()


def fasta_text(a: Analysis) -> str:
    """The contigs as FASTA text (`\\n` line endings). A contig with a letter outside ASCII is
    refused."""
    for contig in a.contigs:
        if not contig.sequence.isascii():
            raise EuleraseError(
                f"Contig ctg{contig.id} holds letters outside plain ASCII, which a FASTA file "
                "cannot hold; change the input letters to ASCII."
            )
    return _written(a, "fasta")


def genbank_text(a: Analysis) -> str:
    """The contigs as GenBank text, with each ORF a CDS feature and the date field 01-JAN-1980.
    Only DNA is written: other letters are refused."""
    if not a.dna:
        raise EuleraseError(
            "GenBank files describe DNA and this input holds other letters; "
            "use the FASTA export, which holds every contig."
        )
    return _written(a, "genbank")


def write_text_file(path: str | os.PathLike, text: str) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)

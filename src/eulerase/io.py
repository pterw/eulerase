"""Reading and writing files: the sequence file reader (FASTA, FASTQ, GenBank, EMBL; plain or
gzip) and the writers (FASTA, GenBank, GFA, CSV)."""

import csv
import gzip
import io
import os
import re
import warnings
import zlib

from Bio import BiopythonParserWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from eulerase.analysis import TABLES, Analysis, columns, data_rows, format_cell, settings_text
from eulerase.assembly import gfa_links
from eulerase.errors import EuleraseError
from eulerase.orfs import Orf, orf_rule_text
from eulerase.sequence import Record, parse_records

_GZIP_MAGIC = b"\x1f\x8b"
_UTF16_MARKS = (b"\xff\xfe", b"\xfe\xff")
_BOM = "\ufeff"
_FIRST_LINE_FORMATS = (("LOCUS ", "genbank"), ("ID ", "embl"), (">", "fasta"), ("@", "fastq"))


def _text_of(data: bytes) -> str:
    """Decode (UTF-16 when the bytes start with its byte-order mark, else UTF-8), drop the
    byte-order mark, read `\\r\\n` and a lone `\\r` as `\\n`, then drop any leading blank lines
    (once, for all users)."""
    encoding = "utf-16" if data.startswith(_UTF16_MARKS) else "utf-8"
    text = data.decode(encoding, errors="replace").removeprefix(_BOM)
    text = text.replace("\r\n", "\n").replace("\r", "\n")
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
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
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


_TOPOLOGY_LINE = (
    "Topology: a cycle in the assembly graph; it can be a circular molecule or a repeat."
)
_SOURCE_NOTE = "assembled by Eulerase from the input; organism not known"


def _cds_location(orf: Orf, length: int) -> SimpleLocation | CompoundLocation:
    """The ORF's location; one that crosses position 1 (`end` above `length`) is two parts in
    transcript order, which on the minus strand is the wrapped part first."""
    if orf.end <= length:
        return SimpleLocation(orf.start, orf.end, strand=orf.strand)
    parts = [
        SimpleLocation(orf.start, length, strand=orf.strand),
        SimpleLocation(0, orf.end - length, strand=orf.strand),
    ]
    return CompoundLocation(parts if orf.strand == 1 else parts[::-1])


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
                f"len={format_cell(contig.length)}",
                f"mean_multiplicity={format_cell(contig.mean_multiplicity)}",
                f"gc={format_cell(contig.gc)} orfs={count}",
                *(["circular=true"] if contig.circular else []),
                settings,
                *rule,
            ]
        )
        name = f"ctg{contig.id}"
        rec = SeqRecord(Seq(contig.sequence), id=name, name=name, description=description)
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["date"] = "01-JAN-1980"
        rec.annotations["topology"] = "circular" if contig.circular else "linear"
        topology = [_TOPOLOGY_LINE] if contig.circular else []
        rec.annotations["comment"] = "\n".join([settings, *rule, *topology])
        rec.features.append(
            SeqFeature(
                SimpleLocation(0, len(contig.sequence)),
                type="source",
                qualifiers={"mol_type": ["other DNA"], "note": [_SOURCE_NOTE]},
            )
        )
        for orf in orf_list or ():
            rec.features.append(
                SeqFeature(
                    _cds_location(orf, len(contig.sequence)),
                    type="CDS",
                    qualifiers={
                        "translation": [orf.protein],
                        "product": ["predicted open reading frame"],
                        "note": [f"frame {orf.frame:+d}"],
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
    """The contigs as FASTA text (`\\n` line endings). A contig with a letter outside ASCII, or
    with `>`, is refused."""
    for contig in a.contigs:
        if not contig.sequence.isascii():
            raise EuleraseError(
                f"Contig ctg{contig.id} holds letters outside plain ASCII, which a FASTA file "
                "cannot hold; change the input letters to ASCII."
            )
        if ">" in contig.sequence:
            raise EuleraseError(
                f"Contig ctg{contig.id} holds '>', which a FASTA file would read back as a new "
                "record; remove '>' from the input letters."
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


_GFA_SEQUENCE = re.compile(r"[A-Za-z=.]+")


def gfa_text(a: Analysis) -> str:
    """The contigs and their overlap links as GFA 1.0 text (`\\n` line endings); the settings
    travel in the header line. A contig with a letter a GFA file cannot hold is refused."""
    for c in a.contigs:
        if not _GFA_SEQUENCE.fullmatch(c.sequence):
            raise EuleraseError(
                f"Contig ctg{c.id} holds characters a GFA file cannot hold (only the letters "
                "A to Z in either case, '=' and '.' are allowed); change the input letters."
            )
    lines = [f"H\tVN:Z:1.0\tst:Z:{settings_text(a.settings)}"]
    for c in a.contigs:
        lines.append(f"S\tctg{c.id}\t{c.sequence}\tLN:i:{c.length}\tKC:i:{c.kmer_count}")
    cycles = {c.id for c in a.contigs if c.circular}
    for src, src_orient, dst, dst_orient in gfa_links(a.contigs, a.k, a.settings.both_strands):
        overlap = 0 if src in cycles else a.k - 1  # a cycle's one link is to itself
        lines.append(f"L\tctg{src}\t{src_orient}\tctg{dst}\t{dst_orient}\t{overlap}M")
    return "\n".join(lines) + "\n"


def csv_text(a: Analysis, table: str) -> str:
    """One table as CSV text: a header row, then every row (none cut), `\\n` after each."""
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(columns(table))
    writer.writerows(data_rows(a, table))  # CSV files hold data only
    return out.getvalue()


def csv_paths(prefix: str) -> dict[str, str]:
    """The five CSV file paths for `prefix`, by lower-case table name, in table order."""
    return {table.lower(): f"{prefix}_{table.lower()}.csv" for table in TABLES}


def write_text_file(path: str | os.PathLike, text: str) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)

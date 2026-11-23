"""Reading and writing files: sequence file reader (FASTA, FASTQ, GenBank, EMBL; plain or gzip)."""

import gzip
import os
import zlib

from Bio import SeqIO

from eulerase.errors import EuleraseError
from eulerase.sequence import Record

MAX_LOAD_CHARS = 60_000

_FORMATS = {
    ".fq": "fastq",
    ".fastq": "fastq",
    ".gb": "genbank",
    ".gbk": "genbank",
    ".genbank": "genbank",
    ".embl": "embl",
}


def read_records(path: str | os.PathLike) -> list[Record]:
    """Read FASTA / FASTQ / GenBank / EMBL (including .gz) with Biopython."""
    name = os.fspath(path)
    if not os.path.exists(name):
        raise EuleraseError(f"The file '{name}' does not exist; check the path.")
    lower = name.lower()
    is_gz = lower.endswith(".gz")
    if is_gz:
        lower = lower[:-3]
    fmt = _FORMATS.get(os.path.splitext(lower)[1], "fasta")
    records: list[Record] = []
    total = 0
    try:
        opener = gzip.open if is_gz else open
        with opener(name, "rt", encoding="utf-8", errors="replace") as handle:
            for rec in SeqIO.parse(handle, fmt):
                s = str(rec.seq).upper()
                if total + len(s) > MAX_LOAD_CHARS:
                    break
                records.append(Record(rec.id, s))
                total += len(s)
    except OSError as err:
        raise EuleraseError(
            f"The file '{name}' could not be read ({err.strerror or 'unreadable'}); "
            "check the path and permissions."
        ) from None
    except (EOFError, zlib.error):
        raise EuleraseError(
            f"The file '{name}' is damaged or incomplete; download or copy it again."
        ) from None
    except ValueError:
        raise EuleraseError(
            f"The file '{name}' could not be read as {fmt.upper()}; check that it is a valid file."
        ) from None
    if not records:
        raise EuleraseError(
            f"No usable {fmt.upper()} records found in '{name}' "
            f"(or the first one exceeds {MAX_LOAD_CHARS:,} bases)."
        )
    return records

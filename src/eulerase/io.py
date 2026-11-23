"""Reading and writing files: sequence file reader (FASTA, FASTQ, GenBank, EMBL; plain or gzip)."""

import gzip
import io
import os
import re
import zlib

from Bio import SeqIO

from eulerase.errors import EuleraseError
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

"""Reading sequence files: every format, plain or gzip, loads every record."""

import codecs
import gzip
import io

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from genomes import random_dna

from eulerase import kmers
from eulerase.errors import EuleraseError
from eulerase.io import detect_format, read_records
from eulerase.sequence import Record

_WRITTEN = [
    Record("alpha", random_dna(90, 1)),
    Record("beta", random_dna(47, 2)),
    Record("gamma", random_dna(130, 3)),
]


def _bio_text(records, fmt):
    out = io.StringIO()
    for r in records:
        rec = SeqRecord(Seq(r.sequence), id=r.name, name=r.name, description="")
        if fmt == "fastq":
            rec.letter_annotations["phred_quality"] = [30] * len(r.sequence)
        elif fmt in ("genbank", "embl"):
            rec.annotations["molecule_type"] = "DNA"
        SeqIO.write(rec, out, fmt)
    return out.getvalue()


def _assert_same_records(got, expected):
    assert [(r.name, r.sequence) for r in got] == [(r.name, r.sequence) for r in expected]


@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gzip"])
@pytest.mark.parametrize(
    ("fmt", "ext"),
    [("fasta", ".fasta"), ("fastq", ".fq"), ("genbank", ".gb"), ("embl", ".embl")],
)
def test_each_format_loads_every_sequence(tmp_path, fmt, ext, gz):
    path = tmp_path / ("reads" + ext + (".gz" if gz else ""))
    text = _bio_text(_WRITTEN, fmt)
    if gz:
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)
    _assert_same_records(read_records(path), _WRITTEN)


def test_oracle_rejects_dropped_sequences():
    emptied = [_WRITTEN[0], Record("beta", ""), _WRITTEN[2]]
    with pytest.raises(AssertionError):
        _assert_same_records(emptied, _WRITTEN)


def test_missing_path_is_an_error_not_a_sequence(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(EuleraseError):
        read_records("ACGT")


def test_a_directory_is_an_error_not_an_oserror(tmp_path):
    with pytest.raises(EuleraseError):
        read_records(tmp_path)


@pytest.mark.filterwarnings("ignore::Bio.BiopythonParserWarning")
def test_biopython_failure_becomes_one_sentence(tmp_path):
    path = tmp_path / "bad.gb"
    path.write_text(
        "LOCUS       X1                10 bp    DNA     linear   UNK 01-JAN-1980\n"
        "FEATURES             Location/Qualifiers\n"
        "     CDS             join(1..\n"
        "ORIGIN\n"
        "        1 acgtacgtac\n//\n"
    )
    with pytest.raises(EuleraseError) as caught:
        read_records(path)
    assert "\n" not in str(caught.value)
    assert "Traceback" not in str(caught.value)


_FA = ">r1\nACGTACGTAC\n>r2\nGGGTTTAAAC\n"
_FA_SEQS = ["ACGTACGTAC", "GGGTTTAAAC"]
_UTF16_LE = codecs.BOM_UTF16_LE + _FA.encode("utf-16-le")
_UTF16_BE = codecs.BOM_UTF16_BE + _FA.encode("utf-16-be")
assert _UTF16_LE != _UTF16_BE


@pytest.mark.parametrize(
    ("name", "data", "expected"),
    [
        pytest.param("a.fa", b"\xef\xbb\xbf" + _FA.encode(), _FA_SEQS, id="bom"),
        pytest.param("a.fa", b"\n\n  \n" + _FA.encode(), _FA_SEQS, id="blank-lines-fasta"),
        pytest.param(
            "a.fq",
            b"\n" + _bio_text([Record("q", "ACGTACGTAC")], "fastq").encode(),
            ["ACGTACGTAC"],
            id="blank-lines-fastq",
        ),
        pytest.param("a.fa", b"\xef\xbb\xbf\n" + _FA.encode(), _FA_SEQS, id="bom-then-blank-line"),
        pytest.param("a.txt", _FA.encode(), _FA_SEQS, id="fasta-named-txt"),
        pytest.param(
            "a.fa",
            _bio_text([Record("q", "ACGTACGTAC")], "fastq").encode(),
            ["ACGTACGTAC"],
            id="fastq-named-fa",
        ),
        pytest.param("a.txt", b"ACGTACGTAC\nGGGTTTAAAC\n", _FA_SEQS, id="bare"),
        pytest.param("a.fa", b"", None, id="empty"),
        pytest.param("a.fa", b"hello, world 42\n", None, id="unrecognised"),
        pytest.param("reads.txt", gzip.compress(_FA.encode()), _FA_SEQS, id="gz-txt"),
        pytest.param("reads.gz", _FA.encode(), _FA_SEQS, id="plain-gz"),
        pytest.param("a.fa", _UTF16_LE, _FA_SEQS, id="utf16-le"),
        pytest.param("a.fa", _UTF16_BE, _FA_SEQS, id="utf16-be"),
        pytest.param("a.fa", gzip.compress(_UTF16_LE), _FA_SEQS, id="utf16-le-gz"),
        pytest.param("a.fa", _UTF16_LE + b"\x00", None, id="utf16-odd-length"),
        pytest.param("a.txt", b"LOCUST\nACGT\n", ["LOCUST", "ACGT"], id="locust-lines"),
    ],
)
def test_format_is_detected_from_content(tmp_path, name, data, expected):
    path = tmp_path / name
    path.write_bytes(data)
    if expected is None:
        with pytest.raises(EuleraseError):
            read_records(path)
    else:
        assert [r.sequence for r in read_records(path)] == expected


@pytest.mark.parametrize(
    ("data", "expected"),
    [
        (b">r1\nACGT\n", "fasta"),
        (b"@q\nACGT\n+\nIIII\n", "fastq"),
        (b"LOCUS       X1   4 bp    DNA\nORIGIN\n", "genbank"),
        (b"ID   X1; SV 1; linear; DNA; STD; UNC; 4 BP.\n", "embl"),
        (b"LOCUST\n", "bare"),
        (b"ACGT\nGGTT\n", "bare"),
        (b"acgt\n", "bare"),
        (b"\xef\xbb\xbf>r1\nACGT\n", "fasta"),
        (b"\n\n  \n>r1\nACGT\n", "fasta"),
        (b">r1\r\nACGT\r\n", "fasta"),
        (b">r1\rACGT\r", "fasta"),
    ],
    ids=[
        "fasta",
        "fastq",
        "genbank",
        "embl",
        "locust-is-bare",
        "bare",
        "bare-lower-case",
        "bom",
        "blank-lines",
        "crlf",
        "cr-only",
    ],
)
def test_detect_format_looks_at_content_only(data, expected):
    assert detect_format(data) == expected


_GZ = gzip.compress(_bio_text(_WRITTEN, "fasta").encode())


@pytest.mark.parametrize(
    "data",
    [_GZ[: len(_GZ) // 2], _GZ[:20] + b"\xff" * 8 + _GZ[28:], b"\x1f\x8b" + b"junk" * 8],
    ids=["truncated", "damaged", "bad-header"],
)
def test_broken_gzip_is_an_error_not_a_crash(tmp_path, data):
    path = tmp_path / "reads.fa.gz"
    path.write_bytes(data)
    with pytest.raises(EuleraseError) as caught:
        read_records(path)
    assert "\n" not in str(caught.value)


def test_case_blank_lines_line_endings_and_wrapping_do_not_change_records(tmp_path):
    clean = tmp_path / "clean.fa"
    clean.write_text(_bio_text(_WRITTEN, "fasta"))
    messy = tmp_path / "messy.fa"
    lines = []
    for r in _WRITTEN:
        s = r.sequence.lower()
        lines += [f">{r.name}", "", *(s[i : i + 17] for i in range(0, len(s), 17)), ""]
    messy.write_bytes("\r\n".join(lines).encode())
    _assert_same_records(read_records(messy), read_records(clean))


_FASTA_LINES = [">r1", "ACGTAC", "GGTTAA", ">r2", "TTGCA"]
_BARE_LINES = ["ACGTAC", "GGTTAA"]


@pytest.mark.parametrize(
    ("lines", "fmt", "expected"),
    [
        (_FASTA_LINES, "fasta", [("r1", "ACGTACGGTTAA"), ("r2", "TTGCA")]),
        (_BARE_LINES, "bare", [("seq1", "ACGTAC"), ("seq2", "GGTTAA")]),
    ],
    ids=["fasta-wrapped", "bare"],
)
@pytest.mark.parametrize(
    "endings",
    [("\n",), ("\r\n",), ("\r",), ("\r", "\r\n", "\n")],
    ids=["lf", "crlf", "cr", "mixed"],
)
def test_every_line_ending_gives_the_same_records(tmp_path, lines, fmt, expected, endings):
    data = "".join(
        line + endings[i % len(endings)] for i, line in enumerate(lines)
    ).encode()  # a trailing ending too
    path = tmp_path / "reads.txt"
    path.write_bytes(data)
    assert detect_format(data) == fmt
    assert [(r.name, r.sequence) for r in read_records(path)] == expected


def test_files_over_60000_bases_load_whole(tmp_path):
    written = [Record("big", random_dna(40_000, 4)), Record("more", random_dna(30_000, 5))]
    path = tmp_path / "long.fa"
    path.write_text(_bio_text(written, "fasta"))
    _assert_same_records(read_records(path), written)


def test_a_file_over_the_limit_is_refused_with_the_limit_named(tmp_path, monkeypatch):
    written = [Record("big", random_dna(40_000, 6)), Record("more", random_dna(30_000, 7))]
    path = tmp_path / "long.fa"
    path.write_text(_bio_text(written, "fasta"))
    monkeypatch.setattr(kmers, "MAX_KMERS", 50_000)
    got = read_records(path)
    _assert_same_records(got, written)
    with pytest.raises(EuleraseError, match="50,000"):
        kmers.from_sequences(got, 3, both_strands=False, unique=False)

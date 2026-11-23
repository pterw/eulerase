"""Writing FASTA and GenBank: each file is read back with Biopython, which shares no code with
the writers, and the worked cases are compared byte for byte with texts typed by hand."""

import io
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from genomes import random_dna

from eulerase.analysis import Analysis, Settings, analyze, settings_text
from eulerase.errors import EuleraseError
from eulerase.io import fasta_text, genbank_text, write_text_file
from eulerase.orfs import MIN_ORF_AA
from eulerase.sequence import Record


def _analysis(sequences: list[str], **settings) -> Analysis:
    records = [Record(f"read{i}", s) for i, s in enumerate(sequences, 1)]
    return analyze(records, Settings(**settings))


def _three_reads_with_orfs() -> Analysis:
    # three unrelated 400-letter reads at k=15 share no 14-mer, so each is its own contig
    reads = [random_dna(400, seed) for seed in (11, 12, 13)]
    a = _analysis(reads, k=15)
    assert len(a.contigs) == 3 and sum(len(o) for o in a.orfs) > 0
    return a


def _repeat() -> Analysis:
    return _analysis(["ATGCATGCATGC"], k=3)  # one circular contig, ATGCAT


def _ten_copies() -> Analysis:
    return _analysis(["ATGC"] * 10, k=3)  # one contig, ATGC


def _words(text: str) -> str:
    return " ".join(text.split())


def _parse(text: str, fmt: str) -> list:
    return list(SeqIO.parse(io.StringIO(text), fmt))


@pytest.mark.parametrize("case", [_three_reads_with_orfs, _repeat], ids=["three-reads", "repeat"])
def test_fasta_and_genbank_round_trip(case):
    a = case()
    for fmt, text in (("fasta", fasta_text(a)), ("genbank", genbank_text(a))):
        parsed = _parse(text, fmt)
        assert [r.id for r in parsed] == [f"ctg{c.id}" for c in a.contigs]
        assert [str(r.seq) for r in parsed] == [c.sequence for c in a.contigs]
    for rec, contig, found in zip(
        _parse(genbank_text(a), "genbank"), a.contigs, a.orfs, strict=True
    ):
        assert rec.annotations["topology"] == ("circular" if contig.circular else "linear")
        cds = [f for f in rec.features if f.type == "CDS"]
        assert [(f.location.start, f.location.end, f.location.strand) for f in cds] == [
            (o.start, o.end, o.strand) for o in found
        ]
        for feature in cds:
            made = feature.extract(rec.seq).translate(to_stop=True)
            assert str(made) == feature.qualifiers["translation"][0]


def test_settings_travel_in_fasta_and_genbank():
    a = _analysis(
        [random_dna(60, 3)], k=5, both_strands=True, unique=True, min_multiplicity=1, clip_tips=True
    )
    line = settings_text(a.settings)
    fasta = _parse(fasta_text(a), "fasta")
    gb = _parse(genbank_text(a), "genbank")
    assert all(line in r.description for r in fasta)
    assert all(line in _words(r.description) for r in gb)
    assert all(line in _words(r.annotations["comment"]) for r in gb)


def test_genbank_date_is_set_by_the_writer():
    locus = genbank_text(_repeat()).splitlines()[0]
    assert locus.startswith("LOCUS") and locus.endswith("01-JAN-1980")


@pytest.mark.parametrize(
    "case", [_three_reads_with_orfs, _ten_copies], ids=["with-orfs", "zero-orfs"]
)
def test_orf_rule_is_written_with_orfs(case):
    a = case()
    assert a.orfs is not None
    line = settings_text(a.settings)
    gb = _parse(genbank_text(a), "genbank")
    descriptions = [r.description for r in _parse(fasta_text(a), "fasta")]
    comments = [_words(r.annotations["comment"]) for r in gb]
    assert len(descriptions) == len(comments) == len(a.contigs)
    for written in descriptions + comments:
        rule = written.split(line, 1)[1]  # what follows the settings: the ORF rule
        assert "ATG" in rule and "standard" in rule and str(MIN_ORF_AA) in rule


def test_non_dna_is_written_as_fasta_without_a_rule_and_refused_as_genbank():
    a = _analysis(["LEMONLEMONADE"], k=3)
    assert a.orfs is None
    fasta = _parse(fasta_text(a), "fasta")
    assert len(fasta) == len(a.contigs) > 0
    assert all("orfs=not available" in r.description for r in fasta)
    assert all(r.description.endswith(settings_text(a.settings)) for r in fasta)
    with pytest.raises(EuleraseError, match="DNA") as refused:
        genbank_text(a)
    assert "FASTA" in str(refused.value)


def test_non_ascii_contig_is_refused_as_fasta_with_its_name():
    a = _analysis(["αβγαβγαβ"], k=3)
    with pytest.raises(EuleraseError, match="ctg1") as refused:
        fasta_text(a)
    assert "ASCII" in str(refused.value)


def _first_difference(expected: str, actual: str) -> int | None:
    """Index of the first differing character, None when the texts are identical."""
    if expected == actual:
        return None
    pairs = zip(expected, actual, strict=False)
    return next((i for i, (x, y) in enumerate(pairs) if x != y), min(len(expected), len(actual)))


# The expected texts below are typed from the frozen formats and worked by hand (spec 5.1).
#
# Ten copies of ATGC at k=3: one contig ATGC, two distinct k-mers (ATG, TGC) seen ten times
# each. length 4; edges 2; k-mer count 20; mean multiplicity (10 + 10) / 2 = 10 -> "10.0000";
# gc = (G + C) / (A + C + G + T) = 2 / 4 = 0.5 -> "0.5000"; ATGC has no stop codon, so no ORF:
# orfs=0 and, because the input is DNA, the ORF rule follows the settings.
#
# The repeat ATGCATGCATGC at k=3: ATG and TGC three times, GCA and CAT twice. One circular
# contig ATGCAT: length 6; edges 4; mean multiplicity (3 + 3 + 2 + 2) / 4 = 2.5 -> "2.5000";
# gc = 2 / 6 = 0.3333; no stop codon in ATGCAT, orfs=0.
#
# FASTA: ">" id, one space, the description, newline, the sequence (under 60 letters, one
# line), newline. GenBank: LOCUS columns per the GenBank flat file: name from column 13,
# length ending in column 40, "bp" in 42-43, molecule type in 48-50, topology in 56-63, the
# division UNK in 65-67, the date in 69-79. DEFINITION and COMMENT lines wrap at 79
# characters on a space, continuation lines indented 12; the writer adds a full stop to the
# DEFINITION. No features, no ORFs. ORIGIN lines are 9 digits wide, lower-case, in tens.
_SETTINGS = "k=3; both strands: no; set of k-mers: no; minimum multiplicity: none; tip clipping: no"
_RULE = (
    "ORFs start at ATG, use the standard genetic code, end at a stop codon "
    "and are at least 10 amino acids long."
)
_TEN_FASTA = f">ctg1 len=4 cov=10.0000x gc=0.5000 orfs=0 {_SETTINGS} {_RULE}\nATGC\n"
_REPEAT_FASTA = f">ctg1 len=6 cov=2.5000x gc=0.3333 orfs=0 {_SETTINGS} {_RULE}\nATGCAT\n"
_TEN_GENBANK = "\n".join(
    [
        "LOCUS       ctg1                       4 bp    DNA     linear   UNK 01-JAN-1980",
        "DEFINITION  len=4 cov=10.0000x gc=0.5000 orfs=0 k=3; both strands: no; set of",
        "            k-mers: no; minimum multiplicity: none; tip clipping: no ORFs start",
        "            at ATG, use the standard genetic code, end at a stop codon and are",
        "            at least 10 amino acids long..",
        "ACCESSION   ctg1",
        "VERSION     ctg1",
        "KEYWORDS    .",
        "SOURCE      .",
        "  ORGANISM  .",
        "            .",
        "COMMENT     k=3; both strands: no; set of k-mers: no; minimum multiplicity:",
        "            none; tip clipping: no",
        "            ORFs start at ATG, use the standard genetic code, end at a stop",
        "            codon and are at least 10 amino acids long.",
        "FEATURES             Location/Qualifiers",
        "ORIGIN",
        "        1 atgc",
        "//",
        "",
    ]
)
_REPEAT_GENBANK = "\n".join(
    [
        "LOCUS       ctg1                       6 bp    DNA     circular UNK 01-JAN-1980",
        "DEFINITION  len=6 cov=2.5000x gc=0.3333 orfs=0 k=3; both strands: no; set of",
        "            k-mers: no; minimum multiplicity: none; tip clipping: no ORFs start",
        "            at ATG, use the standard genetic code, end at a stop codon and are",
        "            at least 10 amino acids long..",
        "ACCESSION   ctg1",
        "VERSION     ctg1",
        "KEYWORDS    .",
        "SOURCE      .",
        "  ORGANISM  .",
        "            .",
        "COMMENT     k=3; both strands: no; set of k-mers: no; minimum multiplicity:",
        "            none; tip clipping: no",
        "            ORFs start at ATG, use the standard genetic code, end at a stop",
        "            codon and are at least 10 amino acids long.",
        "FEATURES             Location/Qualifiers",
        "ORIGIN",
        "        1 atgcat",
        "//",
        "",
    ]
)


@pytest.mark.parametrize(
    ("case", "fasta", "genbank"),
    [(_ten_copies, _TEN_FASTA, _TEN_GENBANK), (_repeat, _REPEAT_FASTA, _REPEAT_GENBANK)],
    ids=["ten-copies", "repeat"],
)
def test_text_outputs_are_exact_for_the_worked_cases(case, fasta, genbank):
    a = case()
    assert _first_difference(fasta, fasta_text(a)) is None
    assert _first_difference(genbank, genbank_text(a)) is None


def test_oracle_rejects_a_changed_byte():
    a = _ten_copies()
    changed = _TEN_FASTA.replace("ATGC\n", "ATGG\n")
    assert _first_difference(changed, fasta_text(a)) == len(changed) - 2
    assert _first_difference(_TEN_GENBANK.replace("atgc", "atgg"), genbank_text(a)) is not None


def _written_bytes(tmp_path: Path) -> list[bytes]:
    a = _repeat()
    out = []
    for name, text in (("a.fasta", fasta_text(a)), ("a.gb", genbank_text(a))):
        write_text_file(tmp_path / name, text)
        out.append((tmp_path / name).read_bytes())
    return out


@pytest.mark.xfail(
    condition=sys.platform == "win32",
    strict=True,
    reason="R2: default text mode writes CRLF on Windows",
)
def test_written_files_use_lf_line_endings(tmp_path: Path):
    for data in _written_bytes(tmp_path):
        assert b"\r" not in data


def test_written_files_have_no_byte_order_mark(tmp_path: Path):
    for data in _written_bytes(tmp_path):
        assert not data.startswith(b"\xef\xbb\xbf")

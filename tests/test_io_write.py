"""Writing FASTA, GenBank, GFA and CSV: FASTA and GenBank are read back with Biopython, GFA with a
small parser written here and CSV with the csv module, none of which share code with the writers;
the worked cases are compared byte for byte with texts typed by hand."""

import csv
import io
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest
from Bio import SeqIO
from genomes import brute_kmers, child_env, random_dna, revcomp

from eulerase import analysis
from eulerase.analysis import (
    TABLES,
    Analysis,
    Settings,
    analyze,
    columns,
    limit_rows,
    rows,
    settings_text,
)
from eulerase.errors import EuleraseError
from eulerase.io import csv_paths, csv_text, fasta_text, genbank_text, gfa_text, write_text_file
from eulerase.orfs import MIN_ORF_AA, orf_rule_text
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
    with pytest.raises(EuleraseError, match="ctg1"):
        gfa_text(a)
    # digits and commas pass the FASTA writer but not the GFA one; other letters are still written
    commas = _analysis(['A,B"A,B"A'], k=3)
    assert fasta_text(commas)
    with pytest.raises(EuleraseError, match="ctg1"):
        gfa_text(commas)
    assert "LEMON" in gfa_text(_analysis(["LEMONLEMONADE"], k=3))


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


# CSV: a header row, then one row per item, comma separated, "\n" after every row, numbers as
# format_cell writes them. Nothing here needs quoting. The Assembly text stops at its "clipped
# k-mers" row (the settings and the ORF rule are text, not data).
#
# Ten copies: 20 edges in ten trails of two steps (ATG then TGC, each seen ten times); nodes AT
# (out 10, in 0), GC (in 10, out 0) and TG (10 and 10); one contig as worked out above.
_TEN_CSV = {
    "Edges": "\n".join(
        ["trail,step,global step,k-mer,from,to,times seen"]
        + [
            f"{t},{s},{2 * t - 2 + s},{kmer},{src},{dst},10"
            for t in range(1, 11)
            for s, kmer, src, dst in ((1, "ATG", "AT", "TG"), (2, "TGC", "TG", "GC"))
        ]
        + [""]
    ),
    "Nodes": (
        "node,in,out,balance,role\nAT,0,10,10,trail start\nGC,10,0,-10,trail end\nTG,10,10,0,\n"
    ),
    "Trails": "\n".join(["trail,edges,sequence"] + [f"{t},2,ATGC" for t in range(1, 11)] + [""]),
    "Contigs": (
        "contig,length,edges,k-mer count,mean multiplicity,min multiplicity,gc,orfs,circular,"
        "has n,strand\n1,4,2,20,10.0000,10,0.5000,0,no,no,as-assembled\n"
    ),
    "Assembly": (
        "measure,value\ncontigs,1\ntotal length,4\nlongest,4\nshortest,4\nmean length,4.0000\n"
        "n50,4\nl50,1\nn90,4\ngc,0.5000\nmean multiplicity,10.0000\ndistinct k-mers,2\n"
        "total k-mers,20\nremoved k-mers,0\nclipped branches,0\nclipped k-mers,0\n"
    ),
}
# The repeat: one trail of ten steps, ATG TGC GCA CAT ATG TGC GCA CAT ATG TGC, seen 3 3 2 2 3 3 2
# 2 3 3 times; nodes AT (out 3, in 2: the start), GC (in 3, out 2: the end), CA and TG balanced.
_REPEAT_CSV = {
    "Edges": (
        "trail,step,global step,k-mer,from,to,times seen\n"
        "1,1,1,ATG,AT,TG,3\n1,2,2,TGC,TG,GC,3\n1,3,3,GCA,GC,CA,2\n1,4,4,CAT,CA,AT,2\n"
        "1,5,5,ATG,AT,TG,3\n1,6,6,TGC,TG,GC,3\n1,7,7,GCA,GC,CA,2\n1,8,8,CAT,CA,AT,2\n"
        "1,9,9,ATG,AT,TG,3\n1,10,10,TGC,TG,GC,3\n"
    ),
    "Nodes": (
        "node,in,out,balance,role\nAT,2,3,1,trail start\nCA,2,2,0,\nGC,3,2,-1,trail end\n"
        "TG,3,3,0,\n"
    ),
    "Trails": "trail,edges,sequence\n1,10,ATGCATGCATGC\n",
    "Contigs": (
        "contig,length,edges,k-mer count,mean multiplicity,min multiplicity,gc,orfs,circular,"
        "has n,strand\n1,6,4,10,2.5000,2,0.3333,0,yes,no,as-assembled\n"
    ),
    "Assembly": (
        "measure,value\ncontigs,1\ntotal length,6\nlongest,6\nshortest,6\nmean length,6.0000\n"
        "n50,6\nl50,1\nn90,6\ngc,0.3333\nmean multiplicity,2.5000\ndistinct k-mers,4\n"
        "total k-mers,10\nremoved k-mers,0\nclipped branches,0\nclipped k-mers,0\n"
    ),
}


@pytest.mark.parametrize(
    ("case", "fasta", "genbank", "csv_texts"),
    [
        (_ten_copies, _TEN_FASTA, _TEN_GENBANK, _TEN_CSV),
        (_repeat, _REPEAT_FASTA, _REPEAT_GENBANK, _REPEAT_CSV),
    ],
    ids=["ten-copies", "repeat"],
)
def test_text_outputs_are_exact_for_the_worked_cases(case, fasta, genbank, csv_texts):
    a = case()
    assert _first_difference(fasta, fasta_text(a)) is None
    assert _first_difference(genbank, genbank_text(a)) is None
    assert list(csv_texts) == list(TABLES)
    assert list(csv_paths("out/x").items()) == [
        ("edges", "out/x_edges.csv"),
        ("nodes", "out/x_nodes.csv"),
        ("trails", "out/x_trails.csv"),
        ("contigs", "out/x_contigs.csv"),
        ("assembly", "out/x_assembly.csv"),
    ]
    for table, text in csv_texts.items():
        assert _first_difference(text, csv_text(a, table)) is None


def test_oracle_rejects_a_changed_byte():
    a = _ten_copies()
    changed = _TEN_FASTA.replace("ATGC\n", "ATGG\n")
    assert _first_difference(changed, fasta_text(a)) == len(changed) - 2
    assert _first_difference(_TEN_GENBANK.replace("atgc", "atgg"), genbank_text(a)) is not None


def _written_bytes(tmp_path: Path) -> list[bytes]:
    a = _repeat()
    texts = {"a.fasta": fasta_text(a), "a.gb": genbank_text(a), "a.gfa": gfa_text(a)}
    texts |= {f"a_{table}.csv": csv_text(a, table) for table in TABLES}
    out = []
    for name, text in texts.items():
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


# --- GFA and CSV -----------------------------------------------------------------------------


def _parse_gfa(text: str) -> tuple[dict, list[tuple[str, str, dict]], list[tuple]]:
    """Header tags, segments (name, sequence, tags) and links (from, o, to, o, overlap)."""
    assert text.endswith("\n") and "\r" not in text

    def tags(fields):
        return {name: value for name, _type, value in (f.split(":", 2) for f in fields)}

    header, segments, links = {}, [], []
    for fields in (line.split("\t") for line in text[:-1].split("\n")):
        if fields[0] == "H":
            header = tags(fields[1:])
        elif fields[0] == "S":
            segments.append((fields[1], fields[2], tags(fields[3:])))
        elif fields[0] == "L":
            links.append(tuple(fields[1:]))
    return header, segments, links


def _flip(orient: str) -> str:
    return "-" if orient == "+" else "+"


def _canonical(link: tuple) -> tuple:
    """A link and its mirror image (b-, a- for a+, b+) are one GFA edge: keep the smaller."""
    a, ao, b, bo = link
    return min(link, (b, _flip(bo), a, _flip(ao)))


def _brute_links(contigs: dict[int, tuple[str, bool]], k: int, both: bool) -> list[tuple]:
    """Every oriented pair of contigs whose end and start share k-1 letters, found by trying all
    pairs. A contig links to itself only when it is circular. Mirror images are one link."""
    oriented = [
        (cid, o, seq if o == "+" else revcomp(seq))
        for cid, (seq, _) in contigs.items()
        for o in ("+-" if both else "+")
    ]
    n = k - 1
    found = [
        (a, ao, b, bo)
        for a, ao, s in oriented
        for b, bo, t in oriented
        if s[-n:] == t[:n] and ((a, ao) != (b, bo) or contigs[a][1])
    ]
    first = {}
    for link in found:
        first.setdefault(_canonical(link), link)
    return list(first.values())


def _unit(k: int, seed: int, both: bool) -> str:
    """A letter string whose cyclic (k-1)-mers are all different (and, with both strands, differ
    from their reverse complements), so that repeating it three times gives one circular contig."""
    n = 2 * k + 6
    for s in range(seed, seed + 500):
        unit = random_dna(n, s)
        ring = unit + unit
        nodes = {ring[i : i + k - 1] for i in range(n)}
        if len(nodes) == n and not (both and nodes & {revcomp(x) for x in nodes}):
            return unit
    raise AssertionError(f"no circular unit for k={k}, both strands {both}")


_GFA_KINDS = ("branching", "circular", "random")


def _gfa_case(i: int) -> tuple[str, bool, int, list[str]]:
    kind, both, k = _GFA_KINDS[i % 3] if i < 20 else "empty", i % 2 == 1, 3 + (2 * i) % 9
    if kind == "branching":  # a shared head and stem, then two different tails
        head, stem, tail1, tail2 = (random_dna(max(12, 3 * k), i * 10 + j) for j in range(4))
        reads = [head + stem + tail1, head + stem + tail2]
    elif kind == "circular":
        reads = [_unit(k, i, both) * 3]
    elif kind == "empty":  # a minimum multiplicity that removes every k-mer: no contigs
        reads = ["ACGTACGT"]
    else:
        reads = [random_dna(5 * k + 10, i * 10 + j) for j in range(2)]
    return kind, both, k, reads


@pytest.mark.parametrize("i", range(21))
def test_gfa_round_trip(i):
    kind, both, k, reads = _gfa_case(i)
    a = _analysis(reads, k=k, both_strands=both, min_multiplicity=9 if kind == "empty" else None)
    header, segments, links = _parse_gfa(gfa_text(a))
    assert header["VN"] == "1.0"
    assert [(name, seq, int(tags["LN"])) for name, seq, tags in segments] == [
        (f"ctg{c.id}", c.sequence, len(c.sequence)) for c in a.contigs
    ]
    assert {overlap for *_, overlap in links} <= {f"{k - 1}M"}
    written = [(int(src[3:]), so, int(dst[3:]), do) for src, so, dst, do, _ in links]
    expected = _brute_links({c.id: (c.sequence, c.circular) for c in a.contigs}, k, both)
    # the writer picks one orientation of each link and the oracle another: compare the classes
    assert [_canonical(link) for link in written] == [_canonical(link) for link in expected]
    # not vacuous: a branching input has a link between two contigs, a circular one a link to itself
    if kind == "branching":
        assert any(link[0] != link[2] for link in expected)
    if kind == "circular":
        assert any(c.circular for c in a.contigs)
        assert any(link[0] == link[2] for link in expected)
    if kind == "empty":  # nothing to write but the header lines
        assert not a.contigs and segments == [] and links == []
        assert gfa_text(a).count("\n") == 1
        assert csv_text(a, "Contigs") == ",".join(columns("Contigs")) + "\n"


def test_a_circular_contig_links_to_itself():
    # ATGCATGCATGC at k=3 is the one circular contig ATGCAT; its end AT is its start AT
    _, segments, links = _parse_gfa(gfa_text(_repeat()))
    assert [name for name, _, _ in segments] == ["ctg1"]
    assert links == [("ctg1", "+", "ctg1", "+", "2M")]


def test_oracle_rejects_a_wrong_gfa_link():
    k = 3  # an overlap of 2 letters
    contigs = {
        1: ("ACGTAC", False),  # ends AC like the start of 2 and 4, but is linear: no self link
        2: ("ACGG", False),
        3: ("GGAT", False),
        4: ("ATAC", False),
        5: ("CCTCC", True),  # circular: its end CC is its start CC
    }
    right = [(1, "+", 2, "+"), (2, "+", 3, "+"), (3, "+", 4, "+"), (4, "+", 1, "+")]
    right += [(4, "+", 2, "+"), (5, "+", 5, "+")]
    assert _brute_links(contigs, k, False) == right
    assert _brute_links(contigs, k, False) != [(1, "+", 2, "+"), (2, "+", 3, "+")]
    # with both strands, ACG (reverse complement CGT) and CGA (reverse complement TCG)
    strands = {1: ("ACG", False), 2: ("CGA", False)}
    assert _brute_links(strands, k, True) == [(1, "+", 1, "-"), (1, "+", 2, "+"), (2, "-", 2, "+")]


@pytest.mark.xfail(strict=True, reason="A8: RC and DP tags")
def test_gfa_segments_carry_ln_and_kc():
    reads = ["ATGCATGCATGC"]
    a = _analysis(reads, k=3)
    (contig,) = a.contigs
    counts = Counter(brute_kmers(reads, 3, False))
    kmers = [contig.sequence[i : i + 3] for i in range(len(contig.sequence) - 2)]
    kc = sum(counts[m] for m in kmers)
    # the numbers a tag could mix up are all different: length, distinct k-mers, sum, rounded mean
    assert [contig.length, len(set(kmers)), kc, round(contig.mean_multiplicity)] == [6, 4, 10, 2]
    ((_, _, tags),) = _parse_gfa(gfa_text(a))[1]
    assert {name: value for name, value in tags.items() if name != "LN"} == {"KC": str(kc)}


def _read_csv(text: str) -> list[list[str]]:
    return list(csv.reader(io.StringIO(text, newline="")))


@pytest.mark.parametrize("table", TABLES)
@pytest.mark.parametrize(
    "sequence", ["AACGTAA", 'A,B"A,B"A'], ids=["five-edges", "commas-and-quotes"]
)
def test_csv_rows_round_trip_for_every_table(sequence, table, monkeypatch):
    monkeypatch.setattr(analysis, "MAX_TABLE_ROWS", 2)
    a = _analysis([sequence], k=3)
    every_row = rows(a, table)
    # the patched limit would cut the window's rows; the export must not
    assert len(limit_rows(every_row)[0]) == min(2, len(every_row))
    if table == "Edges":
        assert len(every_row) > 2
    expected = [columns(table), *every_row]
    if table == "Assembly":
        assert [row[0] for row in every_row[-2:]] == ["settings", "orf rule"]
        expected = expected[:-2]
    assert _read_csv(csv_text(a, table)) == [list(row) for row in expected]


def test_settings_travel_in_the_gfa_header_and_not_in_csv():
    a = _analysis(
        [random_dna(60, 3)], k=5, both_strands=True, unique=True, min_multiplicity=1, clip_tips=True
    )
    assert a.orfs is not None
    line = settings_text(a.settings)
    header, _, _ = _parse_gfa(gfa_text(a))
    assert header["st"] == line
    texts = {table: csv_text(a, table) for table in TABLES}
    assert [t for t, text in texts.items() if line in text or orf_rule_text() in text] == []
    measures = [row[0] for row in _read_csv(texts["Assembly"])]
    assert "settings" not in measures and "orf rule" not in measures
    assert measures[-1] == "clipped k-mers"


_HASH_CHILD = """
import json, sys
from eulerase.analysis import TABLES, Settings, analyze
from eulerase.io import csv_text, fasta_text, genbank_text, gfa_text
from eulerase.sequence import Record

reads = ["ACGTACGT", "AACCAACC", "TTGACATTGACC"]
a = analyze([Record(f"r{i}", r) for i, r in enumerate(reads, 1)], Settings(k=4))
texts = [fasta_text(a), genbank_text(a), gfa_text(a), *(csv_text(a, t) for t in TABLES)]
sys.stdout.buffer.write(json.dumps(texts).encode("ascii"))
"""


def _outputs_under(hash_seeds: list[int], scratch: Path) -> list[bytes]:
    """Start one child per hash seed together; return each child's standard output, in order."""
    children = []
    try:
        for seed in hash_seeds:
            env = child_env()
            env["PYTHONHASHSEED"] = str(seed)
            env["MPLCONFIGDIR"] = str(scratch / f"mpl{seed}")
            children.append(
                subprocess.Popen(
                    [sys.executable, "-c", _HASH_CHILD],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    env=env,
                )
            )
        outputs = []
        for child in children:
            out, err = child.communicate(timeout=120)
            assert child.returncode == 0, err.decode(errors="replace")
            outputs.append(out)
        return outputs
    finally:
        for child in children:
            if child.poll() is None:
                child.kill()
                child.wait()


def test_text_outputs_agree_across_hash_seeds(tmp_path: Path):
    one, other = _outputs_under([1, 4], tmp_path)
    texts = json.loads(one)
    assert len(texts) == 3 + len(TABLES) and all(texts)
    assert one == other

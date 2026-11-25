"""Writing FASTA, GenBank, GFA and CSV: FASTA and GenBank are read back with Biopython, GFA with a
small parser written here and CSV with the csv module, none of which share code with the writers;
the worked cases are compared byte for byte with texts typed by hand."""

import csv
import io
import json
import random
import re
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
    return _analysis(["ATGCATGCATGC"], k=3)  # one circular contig, ATGC


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


_SENSE_CODONS = [
    a + b + c
    for a in "ACGT"
    for b in "ACGT"
    for c in "ACGT"
    if a + b + c not in {"TAA", "TAG", "TGA"}
]


def _gene(codons: int, seed: int) -> str:
    rng = random.Random(seed)
    return "ATG" + "".join(rng.choice(_SENSE_CODONS) for _ in range(codons)) + "TAA"


def _has_no_repeated_8mer(circle: str) -> bool:
    ring = circle + circle[:7]
    return len({ring[i : i + 8] for i in range(len(circle))}) == len(circle)


def _crossing_cycle(strand: int) -> Analysis:
    """A cycle (a circle read three times at k=9) whose least rotation cuts through a gene on the
    given strand, so that the ORF found there has `end` above the contig's length."""
    for seed in range(300):
        gene = _gene(25, seed)
        gene = gene if strand == 1 else revcomp(gene)
        circle = random_dna(30, seed) + gene + random_dna(30, seed + 1000)
        if not _has_no_repeated_8mer(circle):
            continue
        a = _analysis([circle * 3], k=9)
        (contig,) = a.contigs
        if contig.circular and any(o.strand == strand and o.end > contig.length for o in a.orfs[0]):
            return a
    raise AssertionError(f"no circle with a gene across position 1 on strand {strand}")


def _cds_locations(text: str) -> list[str]:
    return re.findall(r"^     CDS +(\S+)$", text, re.MULTILINE)


@pytest.mark.parametrize("strand", [1, -1], ids=["plus", "minus"])
def test_an_orf_across_position_1_is_a_joined_location_on_each_strand(strand):
    a = _crossing_cycle(strand)
    (contig,) = a.contigs
    text = genbank_text(a)
    (rec,) = _parse(text, "genbank")
    cds = [f for f in rec.features if f.type == "CDS"]
    locations = _cds_locations(text)
    assert len(cds) == len(locations) == len(a.orfs[0])
    crossing = 0
    for orf, feature, location in zip(a.orfs[0], cds, locations, strict=True):
        assert "frame" not in feature.qualifiers
        sign = "+" if orf.strand == 1 else "-"
        assert feature.qualifiers["note"] == [f"frame {sign}{abs(orf.frame)}"]
        made = feature.extract(rec.seq).translate(to_stop=True)
        assert str(made) == orf.protein
        if orf.end <= contig.length:
            assert location == (
                f"{orf.start + 1}..{orf.end}"
                if orf.strand == 1
                else f"complement({orf.start + 1}..{orf.end})"
            )
            continue
        crossing += 1
        tail, head = contig.length - orf.start, orf.end - contig.length
        assert tail != head  # the two parts differ in length, so a swap cannot go unseen
        join = f"join({orf.start + 1}..{contig.length},1..{head})"
        assert location == (join if orf.strand == 1 else f"complement({join})")
        assert len(feature.location.parts) == 2
    assert crossing >= 1


def test_a_full_circle_orf_is_written_and_read_back():
    circle = _gene(24, 5)  # ATG, 24 sense codons, TAA: 78 letters
    assert len(circle) == 78 and _has_no_repeated_8mer(circle)
    a = _analysis([circle * 3], k=9)
    (contig,) = a.contigs
    assert contig.circular and contig.length == 78
    (full,) = [o for o in a.orfs[0] if o.length_nt == contig.length and o.strand == 1]
    (rec,) = _parse(genbank_text(a), "genbank")
    features = [f for f in rec.features if f.type == "CDS"]
    assert len(features) == len(a.orfs[0])
    made = {str(f.extract(rec.seq).translate(to_stop=True)) for f in features}
    assert full.protein in made
    assert full.end > contig.length  # the written circle is cut inside the gene: a joined location
    (cds,) = [f for f in features if int(f.location.parts[0].start) == full.start]
    assert [(int(p.start), int(p.end)) for p in cds.location.parts] == [
        (full.start, contig.length),
        (0, full.end - contig.length),
    ]


def test_an_orf_ending_at_the_last_letter_is_a_simple_location():
    a = _analysis([_gene(30, 1)], k=9)  # the read is only the gene: ATG, 30 sense codons, TAA
    (contig,) = a.contigs
    (orf,) = a.orfs[0]
    assert orf.end == contig.length
    text = genbank_text(a)
    assert _cds_locations(text) == [f"{orf.start + 1}..{orf.end}"]
    (rec,) = _parse(text, "genbank")
    (cds,) = [f for f in rec.features if f.type == "CDS"]
    assert str(cds.extract(rec.seq).translate(to_stop=True)) == orf.protein


@pytest.mark.parametrize("case", [_three_reads_with_orfs, _repeat], ids=["three-reads", "repeat"])
def test_every_genbank_record_has_one_source_feature_and_no_organism_or_table_qualifier(case):
    a = case()
    records = _parse(genbank_text(a), "genbank")
    assert len(records) == len(a.contigs)
    for rec in records:
        sources = [f for f in rec.features if f.type == "source"]
        assert len(sources) == 1 and rec.features[0] is sources[0]
        (source,) = sources
        assert (int(source.location.start), int(source.location.end)) == (0, len(rec.seq))
        assert source.qualifiers["mol_type"] == ["other DNA"]
        assert source.qualifiers["note"] == [
            "assembled by Eulerase from the input; organism not known"
        ]
        for feature in rec.features:
            assert not {"organism", "transl_table"} & set(feature.qualifiers)


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
    digits = _analysis(["AC1AC1AC1"], k=3)  # only digits outside A to Z
    assert fasta_text(digits)
    with pytest.raises(EuleraseError, match="ctg1"):
        gfa_text(digits)
    assert "LEMON" in gfa_text(_analysis(["LEMONLEMONADE"], k=3))


def test_a_contig_holding_a_header_character_is_refused_as_fasta_with_its_name():
    a = _analysis(["ACG>TCA"], k=3)
    assert ">" in a.contigs[0].sequence
    assert a.contigs[0].sequence.isascii()  # so the ASCII check is not what refuses it
    with pytest.raises(EuleraseError, match="ctg1") as refused:
        fasta_text(a)
    assert ">" in str(refused.value)


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
# contig ATGC (the circle written once); length 4; edges 4; mean multiplicity
# (3 + 3 + 2 + 2) / 4 = 2.5 -> "2.5000"; gc = 2 / 4 = 0.5 -> "0.5000"; no stop codon in ATGC,
# orfs=0.
#
# FASTA: ">" id, one space, the description, newline, the sequence (under 60 letters, one
# line), newline. GenBank: LOCUS columns per the GenBank flat file: name from column 13,
# length ending in column 40, "bp" in 42-43, molecule type in 48-50, topology in 56-63, the
# division UNK in 65-67, the date in 69-79. DEFINITION and COMMENT lines wrap at 80
# characters (the indent included) on a space, continuation lines indented 12; the writer adds
# a full stop to the DEFINITION. No features, no ORFs. ORIGIN lines are 9 digits wide,
# lower-case, in tens.
_SETTINGS = "k=3; both strands: no; set of k-mers: no; minimum multiplicity: none; tip clipping: no"
_RULE = (
    "ORFs: start codon ATG, standard genetic code, at least 25 amino acids, the longest per "
    "stop codon, ending in a stop codon (none read off an end or across a codon with N); "
    "an ORF is a reading frame, not evidence of a gene."
)
_TEN_FASTA = f">ctg1 len=4 mean_multiplicity=10.0000 gc=0.5000 orfs=0 {_SETTINGS} {_RULE}\nATGC\n"
_REPEAT_FASTA = (
    f">ctg1 len=4 mean_multiplicity=2.5000 gc=0.5000 orfs=0 circular=true {_SETTINGS} {_RULE}\n"
    "ATGC\n"
)
# GFA: the H line with the settings; one S line per contig with LN (length) and KC (the sum of
# the multiplicities of the contig's k-mers); one L line per overlap of k - 1 = 2 letters.
# Ten copies: LN 4, 2 distinct k-mers, KC 10 + 10 = 20, mean 10: all different; no link.
# The repeat: LN 4, 4 distinct k-mers, KC 3 + 3 + 2 + 2 = 10, mean 2.5; the circular contig
# ATGC has the one link to itself, and a link that names a cycle has the overlap 0M.
_TEN_GFA = f"H\tVN:Z:1.0\tst:Z:{_SETTINGS}\nS\tctg1\tATGC\tLN:i:4\tKC:i:20\n"
_REPEAT_GFA = (
    f"H\tVN:Z:1.0\tst:Z:{_SETTINGS}\nS\tctg1\tATGC\tLN:i:4\tKC:i:10\nL\tctg1\t+\tctg1\t+\t0M\n"
)
# GenBank FEATURES: one source feature spanning the record (1..4), its qualifiers from column 22,
# the note wrapped on a space so that no line passes 79 characters (21 + 58).
_SOURCE_FEATURE = [
    "     source          1..4",
    '                     /mol_type="other DNA"',
    '                     /note="assembled by Eulerase from the input; organism not',
    '                     known"',
]
_TEN_GENBANK = "\n".join(
    [
        "LOCUS       ctg1                       4 bp    DNA     linear   UNK 01-JAN-1980",
        "DEFINITION  len=4 mean_multiplicity=10.0000 gc=0.5000 orfs=0 k=3; both strands:",
        "            no; set of k-mers: no; minimum multiplicity: none; tip clipping: no",
        "            ORFs: start codon ATG, standard genetic code, at least 25 amino",
        "            acids, the longest per stop codon, ending in a stop codon (none read",
        "            off an end or across a codon with N); an ORF is a reading frame, not",
        "            evidence of a gene..",
        "ACCESSION   ctg1",
        "VERSION     ctg1",
        "KEYWORDS    .",
        "SOURCE      .",
        "  ORGANISM  .",
        "            .",
        "COMMENT     k=3; both strands: no; set of k-mers: no; minimum multiplicity:",
        "            none; tip clipping: no",
        "            ORFs: start codon ATG, standard genetic code, at least 25 amino",
        "            acids, the longest per stop codon, ending in a stop codon (none read",
        "            off an end or across a codon with N); an ORF is a reading frame, not",
        "            evidence of a gene.",
        "FEATURES             Location/Qualifiers",
        *_SOURCE_FEATURE,
        "ORIGIN",
        "        1 atgc",
        "//",
        "",
    ]
)
_REPEAT_GENBANK = "\n".join(
    [
        "LOCUS       ctg1                       4 bp    DNA     circular UNK 01-JAN-1980",
        "DEFINITION  len=4 mean_multiplicity=2.5000 gc=0.5000 orfs=0 circular=true k=3;",
        "            both strands: no; set of k-mers: no; minimum multiplicity: none; tip",
        "            clipping: no ORFs: start codon ATG, standard genetic code, at least",
        "            25 amino acids, the longest per stop codon, ending in a stop codon",
        "            (none read off an end or across a codon with N); an ORF is a reading",
        "            frame, not evidence of a gene..",
        "ACCESSION   ctg1",
        "VERSION     ctg1",
        "KEYWORDS    .",
        "SOURCE      .",
        "  ORGANISM  .",
        "            .",
        "COMMENT     k=3; both strands: no; set of k-mers: no; minimum multiplicity:",
        "            none; tip clipping: no",
        "            ORFs: start codon ATG, standard genetic code, at least 25 amino",
        "            acids, the longest per stop codon, ending in a stop codon (none read",
        "            off an end or across a codon with N); an ORF is a reading frame, not",
        "            evidence of a gene.",
        "            Topology: a cycle in the assembly graph; it can be a circular",
        "            molecule or a repeat.",
        "FEATURES             Location/Qualifiers",
        *_SOURCE_FEATURE,
        "ORIGIN",
        "        1 atgc",
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
        "contig,length,edges,k-mer count,mean multiplicity,min multiplicity,gc,orfs,cycle,"
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
        "contig,length,edges,k-mer count,mean multiplicity,min multiplicity,gc,orfs,cycle,"
        "has n,strand\n1,4,4,10,2.5000,2,0.5000,0,yes,no,as-assembled\n"
    ),
    "Assembly": (
        "measure,value\ncontigs,1\ntotal length,4\nlongest,4\nshortest,4\nmean length,4.0000\n"
        "n50,4\nl50,1\nn90,4\ngc,0.5000\nmean multiplicity,2.5000\ndistinct k-mers,4\n"
        "total k-mers,10\nremoved k-mers,0\nclipped branches,0\nclipped k-mers,0\n"
    ),
}


@pytest.mark.parametrize(
    ("case", "fasta", "genbank", "gfa", "csv_texts"),
    [
        (_ten_copies, _TEN_FASTA, _TEN_GENBANK, _TEN_GFA, _TEN_CSV),
        (_repeat, _REPEAT_FASTA, _REPEAT_GENBANK, _REPEAT_GFA, _REPEAT_CSV),
    ],
    ids=["ten-copies", "repeat"],
)
def test_text_outputs_are_exact_for_the_worked_cases(case, fasta, genbank, gfa, csv_texts):
    a = case()
    assert _first_difference(fasta, fasta_text(a)) is None
    assert _first_difference(genbank, genbank_text(a)) is None
    assert _first_difference(gfa, gfa_text(a)) is None
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


def _written(tmp_path: Path) -> list[tuple[str, bytes]]:
    """Each text and the bytes of the file written from it. The Greek input reaches only the CSV
    tables: the other writers refuse it."""
    a = _repeat()
    texts = {"a.fasta": fasta_text(a), "a.gb": genbank_text(a), "a.gfa": gfa_text(a)}
    texts |= {f"a_{table}.csv": csv_text(a, table) for table in TABLES}
    greek = _analysis(["αβγαβγαβ"], k=3)
    texts |= {f"g_{table}.csv": csv_text(greek, table) for table in TABLES}
    assert any(not text.isascii() for text in texts.values())
    out = []
    for name, text in texts.items():
        write_text_file(tmp_path / name, text)
        out.append((text, (tmp_path / name).read_bytes()))
    return out


def test_written_files_use_lf_line_endings(tmp_path: Path):
    write_text_file(tmp_path / "ab.txt", "a\nb\n")
    assert (tmp_path / "ab.txt").read_bytes() == b"a\nb\n"
    for _, data in _written(tmp_path):
        assert b"\r" not in data


def test_written_files_are_utf8_without_a_byte_order_mark(tmp_path: Path):
    for text, data in _written(tmp_path):
        assert not data.startswith(b"\xef\xbb\xbf")
        assert data.decode("utf-8") == text


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
    """Every oriented pair of contigs that are not cycles whose end and start share k-1 letters,
    found by trying all pairs; a contig links to itself by this rule like any other. A cycle
    takes part in no such pair and has exactly one link, to itself, (id, "+", id, "+"), with both
    strands too. Mirror images are one link."""
    oriented = [
        (cid, o, seq if o == "+" else revcomp(seq))
        for cid, (seq, cycle) in contigs.items()
        if not cycle
        for o in ("+-" if both else "+")
    ]
    n = k - 1
    found = [(a, ao, b, bo) for a, ao, s in oriented for b, bo, t in oriented if s[-n:] == t[:n]]
    found += [(cid, "+", cid, "+") for cid, (_, cycle) in contigs.items() if cycle]
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
    cycles = {f"ctg{c.id}" for c in a.contigs if c.circular}
    for src, _, dst, _, overlap in links:  # a link that names a cycle has no overlap
        assert overlap == ("0M" if {src, dst} & cycles else f"{k - 1}M")
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
    # ATGCATGCATGC at k=3 is the one circular contig ATGC, which has one link, to itself
    _, segments, links = _parse_gfa(gfa_text(_repeat()))
    assert [name for name, _, _ in segments] == ["ctg1"]
    assert links == [("ctg1", "+", "ctg1", "+", "0M")]


def test_oracle_rejects_a_wrong_gfa_link():
    k = 3  # an overlap of 2 letters
    contigs = {
        1: ("ACGTAC", False),  # ends AC like the start of 2 and 4, and where it begins: self link
        2: ("ACGG", False),
        3: ("GGAT", False),
        4: ("ATAC", False),
        5: ("CCTCC", True),  # a cycle: no overlap link, its one link is to itself
    }
    right = [(1, "+", 1, "+"), (1, "+", 2, "+"), (2, "+", 3, "+"), (3, "+", 4, "+")]
    right += [(4, "+", 1, "+"), (4, "+", 2, "+"), (5, "+", 5, "+")]
    assert _brute_links(contigs, k, False) == right
    assert _brute_links(contigs, k, False) != [(1, "+", 2, "+"), (2, "+", 3, "+")]
    # with both strands, ACG (reverse complement CGT) and CGA (reverse complement TCG)
    strands = {1: ("ACG", False), 2: ("CGA", False)}
    assert _brute_links(strands, k, True) == [(1, "+", 1, "-"), (1, "+", 2, "+"), (2, "-", 2, "+")]


def test_the_gfa_of_a_contig_whose_ends_overlap_links_it_to_itself_with_the_overlap():
    a = _analysis(["AACGTAA", "AAG"], k=3)
    _, segments, links = _parse_gfa(gfa_text(a))
    assert segments == [
        ("ctg1", "AACGTAA", {"LN": "7", "KC": "5"}),
        ("ctg2", "AAG", {"LN": "3", "KC": "1"}),
    ]
    assert not any(c.circular for c in a.contigs)
    assert set(links) == {("ctg1", "+", "ctg1", "+", "2M"), ("ctg1", "+", "ctg2", "+", "2M")}


@pytest.mark.parametrize("read", ["A=A.A=A.", "a=a.a=a."], ids=["upper", "lower"])
def test_the_gfa_refusal_names_what_it_accepts(read):
    assert "ctg1" in gfa_text(_analysis([read], k=3))  # letters, "=" and "." are written
    for bad in ("AC1AC1AC1", 'A,B"A,B"A', "αβγαβγαβ"):
        with pytest.raises(EuleraseError, match="ctg1") as refused:
            gfa_text(_analysis([bad], k=3))
        assert "'='" in str(refused.value) and "'.'" in str(refused.value)


def test_gfa_segments_carry_ln_and_kc():
    reads = ["ACGTT", "ACGTT"]  # one linear contig ACGTT: ACG, CGT and GTT seen twice each
    a = _analysis(reads, k=3)
    (contig,) = a.contigs
    counts = Counter(brute_kmers(reads, 3, False))
    kmers = [contig.sequence[i : i + 3] for i in range(len(contig.sequence) - 2)]
    kc = sum(counts[m] for m in kmers)
    # the numbers a tag could mix up are all different: length, distinct k-mers, sum, rounded mean
    assert [contig.length, len(set(kmers)), kc, round(contig.mean_multiplicity)] == [5, 3, 6, 2]
    assert not contig.circular
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

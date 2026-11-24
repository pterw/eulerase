"""Tests for the analysis result: the counts, trails, contigs and notes `analyze` returns."""

import dataclasses
import re
from collections import Counter

import pytest
from genomes import brute_kmers, random_dna

from eulerase.analysis import Settings, analyze, settings_text
from eulerase.assembly import Metrics
from eulerase.errors import EuleraseError
from eulerase.orfs import Orf
from eulerase.sequence import Record


def _records(*sequences: str) -> list[Record]:
    return [Record(f"seq{i}", s) for i, s in enumerate(sequences, 1)]


def _spectrum_of(items: list[str]) -> dict[int, int]:
    """Independent spectrum: how many distinct k-mers were seen once, twice, ..."""
    return dict(sorted(Counter(Counter(items).values()).items()))


def test_oracle_rejects_a_wrong_spectrum():
    items = ["AAA", "AAA", "AAC"]  # AAA twice, AAC once: one k-mer seen once, one seen twice
    assert _spectrum_of(items) == {1: 1, 2: 1}
    assert _spectrum_of(items) != {1: 1, 2: 2}  # what counting occurrences instead would give


@pytest.mark.parametrize("both_strands", [False, True])
@pytest.mark.parametrize("unique", [False, True])
def test_analyze_matches_independent_counts(both_strands, unique):
    seqs = [random_dna(60, seed=3), random_dna(45, seed=4)]
    a = analyze(_records(*seqs), Settings(k=3, both_strands=both_strands, unique=unique))
    items = brute_kmers(seqs, 3, both_strands)
    seen = Counter(items)
    distinct = set(items)
    steps = [s for trail in a.trails for s in trail]
    walked = sorted(distinct) if unique else items  # the edges of the walk graph
    assert a.total_kmers == (len(distinct) if unique else len(items))
    assert a.distinct_kmers == len(distinct)
    assert dict(a.spectrum) == _spectrum_of(sorted(distinct) if unique else items)
    assert len(steps) == a.total_kmers
    assert len(a.nodes) == len({x for km in distinct for x in (km[:-1], km[1:])})
    assert sorted(s.kmer for s in steps) == (sorted(distinct) if unique else sorted(items))
    assert all(s.multiplicity == (1 if unique else seen[s.kmer]) for s in steps)
    # Node rows: ascending node order; degrees are counts of walked k-mers by their ends.
    assert [n.node for n in a.nodes] == sorted({x for km in distinct for x in (km[:-1], km[1:])})
    for n in a.nodes:
        assert n.in_degree == sum(km[1:] == n.node for km in walked)
        assert n.out_degree == sum(km[:-1] == n.node for km in walked)
        assert n.balance == n.out_degree - n.in_degree
    # Steps: numbered 1..n over all trails, 1.. within a trail, joined end to end.
    assert [s.global_step for s in steps] == list(range(1, len(steps) + 1))
    for i, trail in enumerate(a.trails, 1):
        assert [s.trail for s in trail] == [i] * len(trail)
        assert [s.step for s in trail] == list(range(1, len(trail) + 1))
    assert all((s.source, s.target) == (s.kmer[:-1], s.kmer[1:]) for s in steps)
    assert len({(s.source, s.target, s.key) for s in steps}) == len(steps)


@pytest.mark.parametrize(
    ("reads", "unique", "n_steps", "spelled", "headline", "contigs", "roles"),
    [
        # ATGC ten times at k=3: k-mers ATG and TGC, ten each, so 20 walk edges over the
        # three nodes AT, TG, GC. Every copy is its own trail spelling ATGC; the assembly
        # graph has two edges of multiplicity 10 (mean 10.0) and its one contig is ATGC.
        # All ten trails start at AT and end at GC; TG is only passed through.
        (
            ["ATGC"] * 10,
            False,
            20,
            ("ATGC",) * 10,
            "No single Eulerian path · minimum cover = 10 trails · 20 edges · 3 nodes",
            [("ATGC", False, 10.0)],
            [("AT", "trail start"), ("GC", "trail end"), ("TG", "")],
        ),
        # The same reads as a set: ATG and TGC once each (mean 1.0), one trail ATG -> TGC.
        (
            ["ATGC"] * 10,
            True,
            2,
            ("ATGC",),
            "Eulerian path · 2 edges · 3 nodes",
            [("ATGC", False, 1.0)],
            [("AT", "trail start"), ("GC", "trail end"), ("TG", "")],
        ),
        # ATGCATGCATGC at k=3: ATG and TGC three times, GCA and CAT twice, 10 edges over
        # AT, TG, GC, CA. AT has out 3 and in 2, GC has out 2 and in 3, so it is a path
        # (not a circuit) and one trail spells the input. The assembly graph is the cycle
        # AT -> TG -> GC -> CA -> AT, read from its smallest node AT: ATGCAT, circular. The
        # multiplicities 3, 3, 2, 2 give mean 10 / 4 = 2.5. The one trail starts at AT and
        # ends at GC.
        (
            ["ATGCATGCATGC"],
            False,
            10,
            ("ATGCATGCATGC",),
            "Eulerian path · 10 edges · 4 nodes",
            [("ATGCAT", True, 2.5)],
            [("AT", "trail start"), ("CA", ""), ("GC", "trail end"), ("TG", "")],
        ),
    ],
    ids=["ten-copies", "ten-copies-as-a-set", "repeat-in-one-sequence"],
)
def test_trails_spell_the_input_for_the_worked_cases(
    reads, unique, n_steps, spelled, headline, contigs, roles
):
    a = analyze(_records(*reads), Settings(k=3, unique=unique))
    assert sum(len(t) for t in a.trails) == n_steps
    assert a.spelled == spelled
    assert a.headline == headline
    assert [(c.sequence, c.circular, c.mean_multiplicity) for c in a.contigs] == contigs
    assert [(n.node, n.role) for n in a.nodes] == roles


def test_headline_names_a_circuit_when_the_trail_closes():
    # ATGCAT at k=3: ATG, TGC, GCA, CAT walk AT -> TG -> GC -> CA -> AT, back to the start.
    # The one trail starts and ends at the same node; no other node is a start or an end.
    a = analyze(_records("ATGCAT"), Settings(k=3))
    assert a.headline == "Eulerian circuit · 4 edges · 4 nodes"
    assert sorted(n.role for n in a.nodes) == ["", "", "", "start & end"]


def test_minimum_multiplicity_everything_removed_is_well_formed():
    dna = random_dna(40, seed=8)
    a = analyze(_records(dna), Settings(k=4, min_multiplicity=99))
    assert a.contigs == ()
    assert a.metrics == Metrics(0, 0, 0, 0, None, 0, 0, 0, None, None)
    assert a.orfs == ()
    assert a.removed_kmers == a.distinct_kmers == len(set(brute_kmers([dna], 4, False)))


def test_non_dna_gives_not_available_not_empty():
    a = analyze(_records("LEMONLEMONADE"), Settings(k=3))
    assert a.dna is False
    assert a.orfs is None
    assert len(a.contigs) >= 1
    assert a.metrics.gc is None
    mixed = analyze(_records("ACGTACGTAC", "LEMONLEMON"), Settings(k=3))  # one record not DNA
    assert mixed.dna is False
    assert mixed.orfs is None
    with pytest.raises(EuleraseError):
        analyze(_records("LEMONLEMONADE"), Settings(k=3, both_strands=True))


def test_result_is_immutable():
    a = analyze(_records("ATGCATGCATGC"), Settings(k=3))
    with pytest.raises(dataclasses.FrozenInstanceError):
        a.k = 5  # type: ignore[misc]
    with pytest.raises(dataclasses.FrozenInstanceError):
        a.settings.k = 5  # type: ignore[misc]
    for name in ("trails", "spelled", "nodes", "contigs", "spectrum", "orfs", "notes"):
        assert isinstance(getattr(a, name), tuple), name
    assert all(isinstance(t, tuple) for t in a.trails + a.spectrum + a.orfs)


# Reads whose k-mers have each been counted by hand at k=4: MAIN has 27 distinct 4-mers,
# BRANCH shares its first 16 letters and then leaves MAIN's path with three new 4-mers
# (ATGC, TGCC, GCCA), which dead-end at node CCA. Five copies of MAIN and one of BRANCH
# make the branch's k-mers seen once and every other k-mer five or six times.
_MAIN = "AGGAGTTAAATCGATGTCTCCTTCTGGCTT"
_BRANCH = "AGGAGTTAAATCGATGCCA"


@pytest.mark.parametrize(
    ("reads", "settings", "counts"),
    [
        (["AAAAAAA", "AC", "GT"], Settings(k=4, unique=True), (2, 3)),  # short, duplicates
        ([_MAIN] * 5 + [_BRANCH], Settings(k=4, min_multiplicity=2), (3,)),  # removed by minimum
        ([_MAIN] * 5 + [_BRANCH], Settings(k=4, clip_tips=True), (1, 3)),  # branch, its k-mers
    ],
    ids=["short-and-duplicates", "minimum-multiplicity", "tip-clipping"],
)
def test_notes_state_every_removal(reads, settings, counts):
    a = analyze(_records(*reads), settings)
    text = " ".join(a.notes)
    for n in counts:
        assert re.search(rf"\b{n}\b", text), (n, a.notes)


def test_minimum_multiplicity_is_applied_before_tip_clipping():
    # Same reads as above. The minimum of 2 removes the three branch k-mers (seen once)
    # first, so no dead-end branch is left for tip clipping: removed 3, clipped 0 and 0.
    # Clipping first would clip them (1 branch, 3 k-mers) and leave nothing for the minimum.
    # The spectrum is of the input, with the 3 branch k-mers still at multiplicity 1.
    reads = [_MAIN] * 5 + [_BRANCH]
    a = analyze(_records(*reads), Settings(k=4, min_multiplicity=2, clip_tips=True))
    assert (a.removed_kmers, a.clipped_branches, a.clipped_kmers) == (3, 0, 0)
    assert dict(a.spectrum) == _spectrum_of(brute_kmers(reads, 4, False))


@pytest.mark.parametrize(
    ("settings", "removed", "n_notes"),
    [
        (Settings(k=4), 0, 0),
        (Settings(k=4, min_multiplicity=1), 0, 0),  # 1 is a legal minimum and removes nothing
        (Settings(k=4, both_strands=True), 0, 1),  # the one note: which strand is reported
    ],
    ids=["defaults", "minimum-one", "both-strands"],
)
def test_notes_only_for_what_was_done(settings, removed, n_notes):
    a = analyze(_records(_MAIN), settings)
    assert a.removed_kmers == removed
    assert len(a.notes) == n_notes


@pytest.mark.parametrize(
    ("reads", "settings", "smaller", "unpaired"),
    [
        # AATC, k=3, both strands: k-mers AAT, ATC and, from the reverse strand GATT, GAT, ATT.
        # Nodes AA->AT->TC and GA->AT->TT meet at AT: the merge at AT splits the graph into four
        # one-edge contigs AAT, ATC, ATT, GAT; AAT/ATT and ATC/GAT are twins, so AAT and ATC are
        # reported.
        (["AATC"], Settings(k=3, both_strands=True), 2, 0),
    ],
    ids=["twin-kept"],
)
def test_note_matches_the_strand_labels(reads, settings, smaller, unpaired):
    a = analyze(_records(*reads), settings)
    labels = Counter(c.strand for c in a.contigs)
    assert (labels["smaller-of-pair"], labels["unpaired"]) == (smaller, unpaired)
    # The strand note says which strand of a pair is reported: it is there only when some contig
    # was one of a pair. The only other note here is the tip clipping one.
    strand_notes = len(a.notes) - (a.clipped_branches > 0)
    assert strand_notes == (1 if smaller else 0)


@pytest.mark.parametrize(
    ("settings", "expected"),
    [
        (
            Settings(),
            "k=3; both strands: no; set of k-mers: no; "
            "minimum multiplicity: none; tip clipping: no",
        ),
        (
            Settings(k=5, both_strands=True),
            "k=5; both strands: yes; set of k-mers: no; "
            "minimum multiplicity: none; tip clipping: no",
        ),
        (
            Settings(k=4, unique=True, min_multiplicity=2, clip_tips=True),
            "k=4; both strands: no; set of k-mers: yes; minimum multiplicity: 2; tip clipping: yes",
        ),
        (
            Settings(k=None, input_mode="kmers"),
            "k=from the k-mers; both strands: no; set of k-mers: no; "
            "minimum multiplicity: none; tip clipping: no",
        ),
    ],
)
def test_settings_text_is_stable(settings, expected):
    assert settings_text(settings) == expected


def test_kmer_list_mode_derives_k():
    a = analyze(_records("ACGT", "CGTA"), Settings(k=None, input_mode="kmers"))
    assert a.k == 4
    assert a.total_kmers == 2
    with pytest.raises(EuleraseError):
        Settings(k=5, input_mode="kmers")
    with pytest.raises(ValueError):
        Settings(input_mode="reads")


@pytest.mark.parametrize(
    ("reads", "settings"),
    [([], Settings(k=3)), (["ACGTACGT"], Settings(k=3, min_multiplicity=0))],
    ids=["no-records", "minimum-below-one"],
)
def test_invalid_requests_are_refused(reads, settings):
    with pytest.raises(EuleraseError):
        analyze(_records(*reads), settings)


def test_orfs_are_found_per_contig_in_contig_order():
    # ATG, 25 distinct sense codons (GCT GCC GCA GCG, CTG CTA CTC CTT, GGT GGC GGA GGG, CCT CCC
    # CCA CCG, GTT GTC GTA GTG, ACT ACC ACA ACG, TCT) and TAA: no in-frame stop before the end and
    # no repeated 11-mer, so one forward ORF of 26 amino acids (M, 4 A, 4 L, 4 G, 4 P, 4 V, 4 T,
    # S), 81 nt with the stop, from 0 to 81. Its reverse complement holds no ORF. The second
    # record, shorter than the first, has none.
    gene = (
        "ATG" + "GCTGCCGCAGCG" + "CTGCTACTCCTT" + "GGTGGCGGAGGG" + "CCTCCCCCACCG"
        + "GTTGTCGTAGTG" + "ACTACCACAACG" + "TCT" + "TAA"
    )  # fmt: skip
    protein = "M" + "AAAA" + "LLLL" + "GGGG" + "PPPP" + "VVVV" + "TTTT" + "S"
    other = "ACGTTGCAAGGCTTAGCCATCCGTACCAGT"
    assert len(gene) == 81 and len(protein) == 26
    a = analyze(_records(gene, other), Settings(k=12))
    assert [c.sequence for c in a.contigs] == [gene, other]
    assert a.orfs == ((Orf(1, 1, 0, 81, 81, 26, protein),), ())

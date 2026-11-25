"""The analysis of one input: walk-graph trails, nodes, contigs, metrics, spectrum, ORFs."""

from collections.abc import Sequence
from dataclasses import dataclass, replace

from eulerase import assembly, graph, kmers, orfs
from eulerase.assembly import Contig, Metrics
from eulerase.errors import EuleraseError
from eulerase.orfs import Orf
from eulerase.sequence import Record

DEFAULT_K = 3
MAX_TABLE_ROWS = 2000

_INPUT_MODES = ("sequences", "kmers")
_CYCLE_NOTE = (
    "A contig marked as a cycle is a closed path in the graph: it can be a circular molecule "
    "or a repeat, and its position 1 is its alphabetically first rotation, not a biological start."
)
_GC_AND_N_NOTE = "GC is computed over A, C, G and T; N is left out of both counts."
_BOTH_STRANDS_NOTE = (
    "Both strands: reported one contig per strand pair (the smaller strand in string order "
    "is kept, which may be the reverse complement of the input)."
)


@dataclass(frozen=True)
class Settings:
    k: int | None = DEFAULT_K
    both_strands: bool = False
    unique: bool = False
    min_multiplicity: int | None = None
    clip_tips: bool = False
    input_mode: str = "sequences"

    def __post_init__(self) -> None:
        if self.input_mode not in _INPUT_MODES:
            raise ValueError(f"input_mode must be one of {_INPUT_MODES}, not {self.input_mode!r}")
        if self.input_mode == "kmers" and self.k is not None:
            raise EuleraseError("With a list of k-mers k is taken from them; leave k unset.")


def settings_text(settings: Settings) -> str:
    """One line naming every setting, in a fixed order; written into exports."""
    k = "from the k-mers" if settings.k is None else settings.k
    minimum = "none" if settings.min_multiplicity is None else settings.min_multiplicity
    return (
        f"k={k}; both strands: {_yes_no(settings.both_strands)}; "
        f"set of k-mers: {_yes_no(settings.unique)}; minimum multiplicity: {minimum}; "
        f"tip clipping: {_yes_no(settings.clip_tips)}"
    )


def _yes_no(flag: bool) -> str:
    return "yes" if flag else "no"


@dataclass(frozen=True)
class Step:
    trail: int
    step: int
    global_step: int
    kmer: str
    source: str
    target: str
    key: int
    multiplicity: int


@dataclass(frozen=True)
class NodeRow:
    node: str
    in_degree: int
    out_degree: int
    balance: int
    role: str


@dataclass(frozen=True)
class Analysis:
    settings: Settings
    k: int
    dna: bool
    total_kmers: int
    distinct_kmers: int
    short_records: int
    headline: str
    trails: tuple[tuple[Step, ...], ...]
    spelled: tuple[str, ...]
    nodes: tuple[NodeRow, ...]
    contigs: tuple[Contig, ...]
    metrics: Metrics
    spectrum: tuple[tuple[int, int], ...]
    removed_kmers: int
    clipped_branches: int
    clipped_kmers: int
    orfs: tuple[tuple[Orf, ...], ...] | None
    notes: tuple[str, ...]


NOT_AVAILABLE = "not available"
TABLES = ("Edges", "Nodes", "Trails", "Contigs", "Assembly")


def format_cell(value: object) -> str:
    """The one number format: digits for int, 4 decimals for float, yes/no for bool."""
    if isinstance(value, bool):  # before int: a bool is an int
        return _yes_no(value)
    if value is None:
        return NOT_AVAILABLE
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def _contig_items(a: Analysis) -> list[tuple[Contig, int | None]]:
    counts = [None] * len(a.contigs) if a.orfs is None else [len(o) for o in a.orfs]
    return list(zip(a.contigs, counts, strict=True))


# The Assembly measures that are text, not data: tables of numbers and CSV files leave them out.
_TEXT_MEASURES = ("settings", "orf rule")


def _assembly_items(a: Analysis) -> list[tuple[str, object]]:
    m = a.metrics
    rule = NOT_AVAILABLE if a.orfs is None else orfs.orf_rule_text()
    return [
        ("contigs", m.contigs),
        ("total length", m.total_length),
        ("longest", m.longest),
        ("shortest", m.shortest),
        ("mean length", m.mean_length),
        ("n50", m.n50),
        ("l50", m.l50),
        ("n90", m.n90),
        ("gc", m.gc),
        ("mean multiplicity", m.mean_multiplicity),
        ("distinct k-mers", a.distinct_kmers),
        ("total k-mers", a.total_kmers),
        ("removed k-mers", a.removed_kmers),
        ("clipped branches", a.clipped_branches),
        ("clipped k-mers", a.clipped_kmers),
        ("settings", settings_text(a.settings)),
        ("orf rule", rule),
    ]


# One definition per table: how to list its items from an Analysis, and (header, getter) pairs.
_TABLE_DEFS = {
    "Edges": (
        lambda a: [s for trail in a.trails for s in trail],
        (
            ("trail", lambda s: s.trail),
            ("step", lambda s: s.step),
            ("global step", lambda s: s.global_step),
            ("k-mer", lambda s: s.kmer),
            ("from", lambda s: s.source),
            ("to", lambda s: s.target),
            ("times seen", lambda s: s.multiplicity),
        ),
    ),
    "Nodes": (
        lambda a: a.nodes,
        (
            ("node", lambda n: n.node),
            ("in", lambda n: n.in_degree),
            ("out", lambda n: n.out_degree),
            ("balance", lambda n: n.balance),
            ("role", lambda n: n.role),
        ),
    ),
    "Trails": (
        lambda a: list(enumerate(zip(a.trails, a.spelled, strict=True), 1)),
        (
            ("trail", lambda t: t[0]),
            ("edges", lambda t: len(t[1][0])),
            ("sequence", lambda t: t[1][1]),
        ),
    ),
    "Contigs": (
        _contig_items,
        (
            ("contig", lambda t: t[0].id),
            ("length", lambda t: t[0].length),
            ("edges", lambda t: t[0].edges),
            ("k-mer count", lambda t: t[0].kmer_count),
            ("mean multiplicity", lambda t: t[0].mean_multiplicity),
            ("min multiplicity", lambda t: t[0].min_multiplicity),
            ("gc", lambda t: t[0].gc),
            ("orfs", lambda t: t[1]),
            ("cycle", lambda t: t[0].circular),
            ("has n", lambda t: t[0].has_n),
            ("strand", lambda t: t[0].strand),
        ),
    ),
    "Assembly": (
        _assembly_items,
        (("measure", lambda t: t[0]), ("value", lambda t: t[1])),
    ),
}


def columns(table: str) -> tuple[str, ...]:
    return tuple(header for header, _ in _TABLE_DEFS[table][1])


def rows(a: Analysis, table: str) -> list[tuple[str, ...]]:
    """Every row of a table, every cell already formatted; nothing is cut."""
    items, cells = _TABLE_DEFS[table]
    return [tuple(format_cell(get(item)) for _, get in cells) for item in items(a)]


def data_rows(a: Analysis, table: str) -> list[tuple[str, ...]]:
    """`rows`, without the text measures of the Assembly table."""
    found = rows(a, table)
    if table == "Assembly":
        return [row for row in found if row[0] not in _TEXT_MEASURES]
    return found


def limit_rows(rows: Sequence[tuple[str, ...]]) -> tuple[list[tuple[str, ...]], str | None]:
    """The first MAX_TABLE_ROWS rows for display, and a note when some were cut."""
    if len(rows) <= MAX_TABLE_ROWS:
        return list(rows), None
    return list(rows[:MAX_TABLE_ROWS]), f"Showing {MAX_TABLE_ROWS} of {len(rows)} rows."


def _steps(g, trails: list[graph.Trail], counts: dict[str, int]) -> tuple[tuple[Step, ...], ...]:
    out: list[tuple[Step, ...]] = []
    global_step = 0
    for ti, trail in enumerate(trails, 1):
        row = []
        for si, (u, v, key) in enumerate(trail, 1):
            global_step += 1
            kmer = g.edges[u, v, key]["kmer"]
            row.append(Step(ti, si, global_step, kmer, u, v, key, counts[kmer]))
        out.append(tuple(row))
    return tuple(out)


def _node_rows(g, trails: list[graph.Trail]) -> tuple[NodeRow, ...]:
    starts = {t[0][0] for t in trails}
    ends = {t[-1][1] for t in trails}
    rows = []
    for n in sorted(g):
        if n in starts and n in ends:
            role = "start & end"
        elif n in starts:
            role = "trail start"
        elif n in ends:
            role = "trail end"
        else:
            role = ""
        rows.append(
            NodeRow(n, g.in_degree(n), g.out_degree(n), g.out_degree(n) - g.in_degree(n), role)
        )
    return tuple(rows)


def _headline(g, trails: list[graph.Trail]) -> str:
    edges, nodes = g.number_of_edges(), len(g)
    if len(trails) != 1:
        cover = f"minimum cover = {len(trails)} trails"
        return f"No single Eulerian path · {cover} · {edges} edges · {nodes} nodes"
    closed = trails[0][0][0] == trails[0][-1][1]
    kind = "Eulerian circuit" if closed else "Eulerian path"
    return f"{kind} · {edges} edges · {nodes} nodes"


def _notes(
    found: kmers.Kmers,
    settings: Settings,
    removed: int,
    branches: int,
    clipped: int,
    contigs: Sequence[assembly.Contig],
    gc: float | None,
    letters: Sequence[str],
) -> tuple[str, ...]:
    notes = []
    if found.short_records:
        notes.append(f"{found.short_records} record(s) shorter than k={found.k} were not used.")
    if found.duplicates_dropped:
        notes.append(f"Set of k-mers: {found.duplicates_dropped} repeated occurrence(s) dropped.")
    if removed:
        what = "1 distinct k-mer" if removed == 1 else f"{removed} distinct k-mers"
        notes.append(f"Removed {what} seen fewer than {settings.min_multiplicity} times.")
    if branches:
        notes.append(
            f"Tip clipping removed {branches} dead-end branch(es), {clipped} distinct k-mers."
        )
    if any(c.strand == "smaller-of-pair" for c in contigs):
        notes.append(_BOTH_STRANDS_NOTE)
    if any(c.circular for c in contigs):
        notes.append(_CYCLE_NOTE)
    if gc is not None and any(c.has_n for c in contigs):
        notes.append(_GC_AND_N_NOTE)
    if letters:
        notes.append(
            "GC and ORFs are not available because the input holds letters other than "
            f"A, C, G, T and N: {', '.join(letters)}."
        )
    return tuple(notes)


def _find_orfs(contigs: Sequence[assembly.Contig]) -> tuple[tuple[Orf, ...], ...]:
    """A cycle is searched as a circle, every other contig as a line."""
    return tuple(orfs.find_orfs(c.sequence, circular=c.circular) for c in contigs)


def analyze(records: Sequence[Record], settings: Settings) -> Analysis:
    """Validate, build the k-mers, the walk graph and its trails, then the assembly graph
    (minimum multiplicity, tip clipping), contigs, metrics, spectrum and ORFs."""
    minimum = settings.min_multiplicity
    if minimum is not None and minimum < 1:
        raise EuleraseError("The minimum multiplicity must be at least 1; give a larger number.")
    if settings.input_mode == "kmers":
        found = kmers.from_kmer_list(
            records, both_strands=settings.both_strands, unique=settings.unique
        )
    else:
        found = kmers.from_sequences(
            records, settings.k, both_strands=settings.both_strands, unique=settings.unique
        )
    k = found.k
    counts = kmers.count_kmers(found.items)
    g = graph.build_walk_graph(found.items)
    trails = graph.trail_cover(g)

    ag = assembly.build_assembly_graph(counts)
    removed = branches = clipped = 0
    if minimum is not None:
        ag, removed = assembly.filter_min_multiplicity(ag, minimum)
    if settings.clip_tips:
        ag, branches, clipped = assembly.clip_tips(ag, k)
    contigs = assembly.extract_contigs(ag, k, both_strands=settings.both_strands)
    metrics = assembly.assembly_metrics(contigs)
    # DNA is decided by the non-empty records: an empty record has no letters to be wrong
    letters = sorted({ch for r in records for ch in r.sequence.upper()} - set("ACGTN"))
    dna = not letters
    if not dna:  # GC is not available: no number may sit beside the note that says so
        contigs = tuple(replace(c, gc=None) for c in contigs)
        metrics = replace(metrics, gc=None)
    return Analysis(
        settings=settings,
        k=k,
        dna=dna,
        total_kmers=len(found.items),
        distinct_kmers=len(counts),
        short_records=found.short_records,
        headline=_headline(g, trails),
        trails=_steps(g, trails, counts),
        spelled=tuple(graph.spell(t) for t in trails),
        nodes=_node_rows(g, trails),
        contigs=contigs,
        metrics=metrics,
        spectrum=tuple(assembly.kmer_spectrum(counts).items()),
        removed_kmers=removed,
        clipped_branches=branches,
        clipped_kmers=clipped,
        orfs=_find_orfs(contigs) if dna else None,
        notes=_notes(found, settings, removed, branches, clipped, contigs, metrics.gc, letters),
    )

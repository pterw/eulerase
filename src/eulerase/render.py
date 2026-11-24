"""Layout and figure of the walk graph: node coordinates and a matplotlib figure, no window."""

import math
import os
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

import networkx as nx
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.text import Text

from eulerase.analysis import Analysis, Step, settings_text
from eulerase.errors import EuleraseError

# ---------------------------------------------------------------------------
# Palette (lemon on neutral surfaces) and limits
# ---------------------------------------------------------------------------
INK = "#1F1F1E"  # text and outlines
MUTED = "#5F5E58"  # secondary text
RULE = "#D9D8D2"  # hairlines and panel borders
SURFACE = "#FCFCFB"  # the canvas and the tables
WINDOW = "#F3F2EE"  # the window's background
LEMON = "#F9D923"  # the primary action and the walk's current position, nothing else
LEMON_SOFT = "#FDF0A6"  # the selected table row
FIELD_BORDER = "#85837A"
PRIMARY_BORDER = "#8C7400"
ERROR = "#B3261E"
FUTURE = "#949289"  # arrows not yet used
TRAIL_COLORS = ("#2A78D6", "#EB6834", "#1BAF7A")
# Trail t takes the colour TRAIL_COLORS[(t - 1) % 3] and the dash TRAIL_DASHES[((t - 1) // 3) % 3],
# so no two of the first nine trails share both.
TRAIL_DASHES = ("solid", "dashed", "dotted")

# One drawing limit with two counts: above either, the graph is not drawn (tables still work).
MAX_DRAW_NODES = 400
MAX_DRAW_EDGES = 2000  # arrows: the walk graph's edges, self-loops and parallel edges included
OVERVIEW_ABOVE_NODES = 60
MAX_COPY_MARKS = 6
FIGURE_ASPECT = 1.6  # width over height of the plot area of a saved figure
FIGURE_FORMATS = ("png", "svg", "pdf")


def _walk_graph(a: Analysis) -> nx.MultiDiGraph:
    """The walk graph again, from the nodes and steps the analysis holds."""
    g = nx.MultiDiGraph()
    g.add_nodes_from(row.node for row in a.nodes)
    for trail in a.trails:
        for s in trail:
            g.add_edge(s.source, s.target, key=s.key, kmer=s.kmer)
    return g


def _counts_over_limit(g: nx.MultiDiGraph) -> list[str]:
    """Each count of the drawing limit that `g` is over, with its number and its limit."""
    over = []
    if len(g) > MAX_DRAW_NODES:
        over.append(f"{len(g)} nodes, over the limit of {MAX_DRAW_NODES}")
    if g.number_of_edges() > MAX_DRAW_EDGES:
        over.append(f"{g.number_of_edges()} arrows, over the limit of {MAX_DRAW_EDGES}")
    return over


def compute_layout(
    a: Analysis, aspect: float = FIGURE_ASPECT
) -> dict[str, tuple[float, float]] | None:
    """Reading-order layout: nodes in the order the walk first meets them, on a back-and-forth
    grid (left to right on one row, right to left on the next) filling a plot area `aspect` wide
    and 1 high.

    The order is: trails in order; in each, the source of its first step, then the target of
    every step; a node already met is skipped. With n nodes, c columns and r = ceil(n / c) rows
    the spacing is min(aspect / c, 1 / r); the grid has the c from 1 to n with the widest spacing,
    the smaller c on a tie. Node i is in row i // c, in column i % c on an even row and
    c - 1 - i % c on an odd one, at x = (column + 0.5) * aspect / c and y = 1 - (row + 0.5) / r.
    There is no iteration and no random number, so the coordinates depend only on the input and
    `aspect`. Returns None above the drawing limit and {} for no nodes; an `aspect` that is not
    finite or is at or below 0 raises ValueError.
    """
    if not (math.isfinite(aspect) and aspect > 0):
        raise ValueError(f"aspect must be a finite number above 0, not {aspect!r}")
    if _counts_over_limit(_walk_graph(a)):
        return None
    order: dict[str, None] = {}
    for trail in a.trails:
        order.setdefault(trail[0].source)
        for step in trail:
            order.setdefault(step.target)
    n = len(order)
    if n == 0:
        return {}
    c = max(range(1, n + 1), key=lambda cols: min(aspect / cols, 1 / math.ceil(n / cols)))
    r = math.ceil(n / c)
    layout = {}
    for i, node in enumerate(order):
        row, offset = divmod(i, c)
        column = offset if row % 2 == 0 else c - 1 - offset
        layout[node] = ((column + 0.5) * aspect / c, 1 - (row + 0.5) / r)
    return layout


# ---------------------------------------------------------------------------
# Label placement: by character counts in a monospaced font, never by font metrics, so the
# places are the same on every platform. Lengths are in points; a box is (x0, y0, x1, y1).
# ---------------------------------------------------------------------------
MONO = "DejaVu Sans Mono"  # ships with matplotlib
MONO_ADVANCE = 0.602  # width of one character of MONO, in em
LINE_HEIGHT = 1.25  # height of a line of MONO, in em, with room for the highest and lowest glyph
LABEL_PAD = 3.0  # room kept around every label, beyond its text
LABEL_PT, WORD_PT = 8.5, 7.5  # font sizes of a k-mer label and of a start or end word
MARKS_LINE_H = 7.0  # height of the line of copy marks under a label of a k-mer seen more than once
Box = tuple[float, float, float, float]
Point = tuple[float, float]


class Spot(NamedTuple):
    """A place for a label: its box goes beside the anchor (x, y), on the side of the unit normal
    (nx, ny) and `gap` away from the anchor."""

    x: float
    y: float
    nx: float
    ny: float
    gap: float


@dataclass(frozen=True)
class Wanted:
    """A label to place: its text, the spots to try in order, the size of one character and of one
    line, and whether it has a line of copy marks under its text. Its full box holds that line;
    a label that fits only without it is shortened. `kind` is "label" (the k-mer of an arrow) or
    "word" (the "start 1" of a ring)."""

    key: str
    text: str
    spots: tuple[Spot, ...]
    char_w: float
    line_h: float
    marks: bool = False
    kind: str = "label"


@dataclass(frozen=True)
class Placed:
    """Where a label went. `box` is what it takes, so that what is drawn under it (the copy marks)
    can be put inside it."""

    text: str
    x: float  # the centre of the text line
    y: float
    box: Box  # what it takes, padding and the marks line included
    marks: bool  # the box holds the marks line, under the text
    kind: str  # of its `Wanted`


@dataclass(frozen=True)
class Placement:
    placed: dict[str, Placed]
    shortened: int  # labels placed without the marks line they wanted
    left_out: int  # labels that fit nowhere


def _crosses(a: Box, b: Box) -> bool:
    return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]


class _Occupied:
    """The boxes already taken, in a grid of cells, so that a check meets only its neighbours."""

    CELL = 48.0

    def __init__(self, boxes: list[Box]) -> None:
        self._cells: dict[tuple[int, int], list[Box]] = defaultdict(list)
        for box in boxes:
            self.add(box)

    def _keys(self, box: Box):
        x0, y0, x1, y1 = (int(v // self.CELL) for v in box)
        return ((i, j) for i in range(x0, x1 + 1) for j in range(y0, y1 + 1))

    def add(self, box: Box) -> None:
        for key in self._keys(box):
            self._cells[key].append(box)

    def hits(self, box: Box) -> bool:
        return any(_crosses(box, o) for key in self._keys(box) for o in self._cells.get(key, ()))


def _box_at(w: Wanted, marks: bool, spot: Spot) -> tuple[float, float, Box]:
    """The centre of the text line and the box of the label at `spot`: beside the anchor, `gap`
    clear of it. The marks line, when kept, lies under the text and inside the box."""
    extra = MARKS_LINE_H if marks else 0.0
    half_w = len(w.text) * w.char_w / 2 + LABEL_PAD
    half_h = (w.line_h + extra) / 2 + LABEL_PAD
    reach = abs(spot.nx) * half_w + abs(spot.ny) * half_h + spot.gap
    x, y = spot.x + spot.nx * reach, spot.y + spot.ny * reach
    return x, y + extra / 2, (x - half_w, y - half_h, x + half_w, y + half_h)


def _first_free(w: Wanted, marks: bool, occupied: _Occupied, bounds: Box):
    """The text centre and box at the first spot of `w` where the label lies inside `bounds` and
    meets no taken box, or None."""
    for spot in w.spots:
        x, y, box = _box_at(w, marks, spot)
        inside = bounds[0] <= box[0] and bounds[1] <= box[1] and box[2] <= bounds[2]
        if inside and box[3] <= bounds[3] and not occupied.hits(box):
            return x, y, box
    return None


def place_labels(wanted: list[Wanted], blocked: list[Box], bounds: Box) -> Placement:
    """Give each label, in order, the first free spot for its full box, else for its text alone
    (shortened), else leave it out. A spot is free when its box lies inside `bounds` and meets
    neither a `blocked` box nor a label already placed. Only the arguments decide the result, so
    a step of the walk never moves a label. Linear in the labels: a check meets only the boxes
    in the cells it covers."""
    occupied = _Occupied(blocked)
    placed: dict[str, Placed] = {}
    cut = 0
    for w in wanted:
        for marks in (True, False) if w.marks else (False,):
            found = _first_free(w, marks, occupied, bounds)
            if found:
                x, y, box = found
                occupied.add(box)
                placed[w.key] = Placed(w.text, x, y, box, marks, w.kind)
                cut += w.marks and not marks
                break
    return Placement(placed, cut, len(wanted) - len(placed))


# ---------------------------------------------------------------------------
# The plot: built once, restyled for a step
# ---------------------------------------------------------------------------
SPOT_GAPS = (3.0, 8.0, 14.0, 21.0)  # how far from its arrow a label may sit, nearest first
SPOT_TS = (0.5, 0.38, 0.62, 0.26, 0.74)  # where along its arrow, the middle first
WORD_GAP = 2.0  # between a start or end word and the outermost ring of its node
BULGES = (0.0, 0.3, -0.3, 0.5, -0.5, 0.75, -0.75)  # bends tried, in node spacings, for an arrow
BEND_CLEARANCE = 3.0  # points an arrow keeps from a node it passes
RING_STEP = 0.22  # each ring at a node lies this fraction of the radius beyond the one before
LOOP_HALF_ANGLE = 0.55  # a self-loop leaves and meets its node this far either side of its middle
LOOP_RAD = 2.0  # how far a self-loop bulges, as the `rad` of an arc
SIBLING_RAD = 0.52  # the bend between two arrows between the same two nodes
REVERSE_RAD = 0.22  # the bend of an arrow whose pair of nodes has an arrow back
NODE_FONTS = (10.0, 8.5, 7.5)  # node letters, when up to 3, up to 5, and more letters fit
OVERVIEW_WEIGHT = 0.5  # line width and arrowhead of an overview arrow, against a full one


@dataclass(frozen=True)
class _Sizes:
    scale: float  # points per unit of the layout
    spacing: float  # between neighbouring nodes, in points
    radius: float  # of a node
    font: float  # of its letters


@dataclass(frozen=True)
class _Plan:
    """One arrow before it is drawn: its steps and its ends and bend in points."""

    name: str  # what the arrow and its label share in their gids
    steps: tuple[tuple[int, int], ...]  # (global step, trail) of each copy
    p0: Point
    p1: Point
    rad: float
    kmer: str | None  # None in the overview, which has no text

    @property
    def gid(self) -> str:
        return f"arrow:{self.name}"

    @property
    def label_gid(self) -> str:
        return f"label:{self.name}"


@dataclass
class _Arrow:
    patch: FancyArrowPatch
    label: Text | None
    steps: tuple[tuple[int, int], ...]
    weight: float


@dataclass
class _Plot:
    arrows: list[_Arrow]
    placement: Placement  # of the labels and words: what WQ draws the marks by, and the counts

    def style(self, upto: int) -> None:
        """Restyle the arrows and their labels for the first `upto` steps of the walk."""
        for arrow in self.arrows:
            used = [trail for step, trail in arrow.steps if step <= upto]
            active = any(step == upto for step, _ in arrow.steps)
            patch = arrow.patch
            patch.set_color(TRAIL_COLORS[(used[0] - 1) % len(TRAIL_COLORS)] if used else FUTURE)
            patch.set_linewidth(arrow.weight * (2.4 if active else 1.6 if used else 1.0))
            patch.set_mutation_scale(arrow.weight * (11 if used else 8))
            patch.set_zorder(3 if used else 1)
            if arrow.label is not None:
                arrow.label.set_color(INK if used else MUTED)
                arrow.label.set_fontweight("bold" if active else "normal")


def _arc_point(p0: Point, p1: Point, rad: float, t: float) -> Point:
    """The point at t on the curve `arc3,rad=rad` from p0 to p1: a quadratic Bezier whose control
    point is the middle of the chord moved by `rad` times the chord turned a quarter."""
    cx = (p0[0] + p1[0]) / 2 + rad * (p1[1] - p0[1])
    cy = (p0[1] + p1[1]) / 2 - rad * (p1[0] - p0[0])
    s = 1 - t
    return (
        s * s * p0[0] + 2 * s * t * cx + t * t * p1[0],
        s * s * p0[1] + 2 * s * t * cy + t * t * p1[1],
    )


def _arrow_spots(p0: Point, p1: Point, rad: float) -> tuple[Spot, ...]:
    """Places for a label along an arrow, nearest first, on the side it bulges to first."""
    d = math.dist(p0, p1) or 1.0
    nx, ny = (p1[1] - p0[1]) / d, -(p1[0] - p0[0]) / d
    first = -1.0 if rad < 0 else 1.0
    return tuple(
        Spot(*_arc_point(p0, p1, rad, t), side * nx, side * ny, gap)
        for gap in SPOT_GAPS
        for t in SPOT_TS
        for side in (first, -first)
    )


def _sizes(layout: dict[str, Point], scale: float, overview: bool) -> _Sizes:
    """The size of a node in points: it holds its letters when the spacing allows and is at most
    a quarter of the spacing, so that neighbours do not touch and a self-loop has room."""
    points = list(layout.values())
    spacing = min((math.dist(p, q) for p, q in zip(points, points[1:])), default=1.0) * scale
    if overview:
        return _Sizes(scale, spacing, max(1.5, min(3.2, 0.2 * spacing)), 0.0)
    letters = len(next(iter(layout)))
    base = NODE_FONTS[0 if letters <= 3 else 1 if letters <= 5 else 2]
    radius = min(0.301 * base * letters + 3.5, 0.25 * spacing)
    font = max(4.0, min(base, (radius - 3.0) / (0.301 * letters)))
    return _Sizes(scale, spacing, radius, font)


def _rings(a: Analysis) -> list[tuple[str, int, str]]:
    """(kind, trail, node) of every ring, trail by trail, the start before the end."""
    return [
        (kind, t, node)
        for t, trail in enumerate(a.trails, 1)
        for kind, node in (("start", trail[0].source), ("end", trail[-1].target))
    ]


def _groups(a: Analysis, overview: bool) -> dict[tuple, list[Step]]:
    """The steps behind each arrow, in the order of the walk: an arrow for each trail and k-mer,
    or in the overview for each pair of nodes."""
    groups: dict[tuple, list[Step]] = {}
    for trail in a.trails:
        for s in trail:
            key = (s.source, s.target) if overview else (s.trail, s.kmer)
            groups.setdefault(key, []).append(s)
    return groups


def _bend(
    p0: Point,
    p1: Point,
    nodes: dict[str, tuple[Point, float]],
    ends: tuple[str, str],
    spacing: float,
) -> float:
    """The first bend of BULGES at which the arrow from p0 to p1 passes clear of every node but
    its two `ends` (each with its centre and the room it takes), or 0 when none does."""
    others = [circle for n, circle in nodes.items() if n not in ends]
    chord = math.dist(p0, p1) or 1.0
    for bulge in BULGES:
        rad = bulge * spacing / chord
        samples = [_arc_point(p0, p1, rad, i / 8) for i in range(1, 8)]
        if all(math.dist(p, q) > room + BEND_CLEARANCE for p in samples for q, room in others):
            return rad
    return 0.0


def _chord_ends(pu: Point, pv: Point, gap_u: float, gap_v: float) -> tuple[Point, Point]:
    d = math.dist(pu, pv) or 1.0
    ux, uy = (pv[0] - pu[0]) / d, (pv[1] - pu[1]) / d
    return (pu[0] + ux * gap_u, pu[1] + uy * gap_u), (pv[0] - ux * gap_v, pv[1] - uy * gap_v)


def _loop_ends(centre: Point, reach: float, index: int) -> tuple[Point, Point]:
    """Where the `index`-th self-loop at a node leaves and meets it: the first goes up, the next
    ones turn a third of a circle each."""
    angle = math.pi / 2 - index * math.pi / 3
    leaves, meets = (angle - LOOP_HALF_ANGLE, angle + LOOP_HALF_ANGLE)
    return (
        (centre[0] + reach * math.cos(leaves), centre[1] + reach * math.sin(leaves)),
        (centre[0] + reach * math.cos(meets), centre[1] + reach * math.sin(meets)),
    )


def _plan_arrows(
    a: Analysis, pos: dict[str, Point], room: dict[str, float], sizes: _Sizes, overview: bool
) -> list[_Plan]:
    """Every arrow's ends and bend, in the order of the walk. Arrows between the same two nodes
    are spread apart, an arrow with one back is bent aside, and a lone straight arrow that would
    pass through another node is bent clear of it."""
    groups = _groups(a, overview)
    between: dict[tuple[str, str], list[tuple]] = defaultdict(list)
    for key, steps in groups.items():
        between[(steps[0].source, steps[0].target)].append(key)
    loops: dict[str, int] = defaultdict(int)
    nodes = {n: (pos[n], room[n]) for n in pos}
    plans = []
    for key, steps in groups.items():
        first = steps[0]
        u, v = first.source, first.target
        if u == v:
            p0, p1 = _loop_ends(pos[u], room[u] + 1.0, loops[u])
            loops[u] += 1
            rad = LOOP_RAD
        else:
            siblings = between[(u, v)]
            rad = SIBLING_RAD * (siblings.index(key) - (len(siblings) - 1) / 2)
            rad += REVERSE_RAD if (v, u) in between else 0.0
            p0, p1 = _chord_ends(pos[u], pos[v], room[u] + 1.5, room[v] + 2.5)
            if rad == 0.0 and not overview:
                rad = _bend(p0, p1, nodes, (u, v), sizes.spacing)
        name = f"{u}:{v}" if overview else f"{first.trail}:{first.kmer}"
        steps_of = tuple((s.global_step, s.trail) for s in steps)
        plans.append(_Plan(name, steps_of, p0, p1, rad, None if overview else first.kmer))
    return plans


def _word(kind: str, trail: int, centre: Point, room: float) -> Wanted:
    """The word of a ring ("start 2"): below its node, else above, else to the right or the left."""
    spots = tuple(
        Spot(centre[0] + nx * room, centre[1] + ny * room, nx, ny, WORD_GAP)
        for nx, ny in ((0.0, -1.0), (0.0, 1.0), (1.0, 0.0), (-1.0, 0.0))
    )
    return Wanted(
        f"word:{kind}:{trail}",
        f"{kind} {trail}",
        spots,
        MONO_ADVANCE * WORD_PT,
        LINE_HEIGHT * WORD_PT,
        kind="word",
    )


def _label(plan: _Plan) -> Wanted:
    """The label of an arrow: its k-mer, along the arrow, with a marks line for several copies."""
    assert plan.kmer is not None
    return Wanted(
        plan.label_gid,
        plan.kmer,
        _arrow_spots(plan.p0, plan.p1, plan.rad),
        MONO_ADVANCE * LABEL_PT,
        LINE_HEIGHT * LABEL_PT,
        len(plan.steps) > 1,
    )


def _add_node(ax, node: str, centre: Point, sizes: _Sizes, overview: bool) -> None:
    x, y = centre[0] / sizes.scale, centre[1] / sizes.scale
    ax.add_patch(
        Circle(
            (x, y),
            sizes.radius / sizes.scale,
            facecolor=INK if overview else SURFACE,
            edgecolor=INK,
            linewidth=0 if overview else 1.5,
            zorder=4,
            gid=f"node:{node}",
        )
    )
    if not overview:
        ax.text(
            x,
            y,
            node,
            ha="center",
            va="center",
            family=MONO,
            fontsize=sizes.font,
            fontweight="bold",
            color=INK,
            zorder=6,
            gid=f"letters:{node}",
        )


def _add_ring(ax, kind: str, trail: int, centre: Point, radius: float, scale: float) -> None:
    ax.add_patch(
        Circle(
            (centre[0] / scale, centre[1] / scale),
            radius / scale,
            fill=False,
            edgecolor=INK,
            linewidth=1.2,
            linestyle=TRAIL_DASHES[0 if kind == "start" else 1],
            zorder=3,
            gid=f"ring:{kind}:{trail}",
        )
    )


def _add_text(ax, key: str, placed: Placed, scale: float) -> Text:
    is_word = placed.kind == "word"
    return ax.text(
        placed.x / scale,
        placed.y / scale,
        placed.text,
        ha="center",
        va="center",
        family=MONO,
        fontsize=WORD_PT if is_word else LABEL_PT,
        color=INK,
        zorder=5,
        gid=key,
        bbox=None
        if is_word
        else dict(boxstyle="round,pad=0.15", facecolor=SURFACE, edgecolor="none", alpha=0.85),
    )


def _add_arrow(ax, plan: _Plan, scale: float, label: Text | None, weight: float) -> _Arrow:
    patch = FancyArrowPatch(
        (plan.p0[0] / scale, plan.p0[1] / scale),
        (plan.p1[0] / scale, plan.p1[1] / scale),
        connectionstyle=f"arc3,rad={plan.rad}",
        arrowstyle="-|>",
        shrinkA=0,
        shrinkB=0,
        gid=plan.gid,
    )
    ax.add_patch(patch)
    return _Arrow(patch, label, plan.steps, weight)


def _build_plot(ax, a: Analysis, layout: dict[str, Point], aspect: float) -> _Plot:
    """Draw the walk graph on `ax`, whose data area is `aspect` wide and 1 high, with no step
    drawn as used yet. Everything a step can restyle exists from the start; `_Plot.style` sets it.
    Above OVERVIEW_ABOVE_NODES nodes the plot is an overview: dots and one arrow for each pair of
    nodes, and no text."""
    overview = len(layout) > OVERVIEW_ABOVE_NODES
    ax.set_xlim(0, aspect)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.apply_aspect()  # the box the axes will have, so its size in points is known now
    if not layout:
        return _Plot([], Placement({}, 0, 0))
    scale = ax.get_position().height * ax.figure.get_figheight() * 72
    sizes = _sizes(layout, scale, overview)
    pos = {node: (x * scale, y * scale) for node, (x, y) in layout.items()}

    rings = [] if overview else _rings(a)
    room = dict.fromkeys(pos, sizes.radius)  # of a node and its rings: nothing may touch it
    ring_radii = []
    for _, _, node in rings:
        room[node] += sizes.radius * RING_STEP
        ring_radii.append(room[node])
    plans = _plan_arrows(a, pos, room, sizes, overview)

    # The labels of the k-mers go first, then the words, in a fixed order: one rule and one count.
    bounds = (0.0, 0.0, aspect * scale, scale)
    blocked = [(x - room[n], y - room[n], x + room[n], y + room[n]) for n, (x, y) in pos.items()]
    wanted = [_label(plan) for plan in plans if plan.kmer is not None]
    wanted += [_word(kind, t, pos[node], room[node]) for kind, t, node in rings]
    placement = place_labels(wanted, blocked, bounds)

    for node, centre in pos.items():
        _add_node(ax, node, centre, sizes, overview)
    for (kind, t, node), radius in zip(rings, ring_radii):
        _add_ring(ax, kind, t, pos[node], radius, scale)
    texts = {key: _add_text(ax, key, text, scale) for key, text in placement.placed.items()}
    weight = OVERVIEW_WEIGHT if overview else 1.0
    arrows = [_add_arrow(ax, plan, scale, texts.get(plan.label_gid), weight) for plan in plans]
    return _Plot(arrows, placement)


def _legend_and_title(fig: Figure, a: Analysis) -> None:
    elements = []
    if len(a.trails) > 1:
        for i, t in enumerate(a.trails, 1):
            c = TRAIL_COLORS[(i - 1) % len(TRAIL_COLORS)]
            elements.append(Line2D([0], [0], color=c, lw=2.2, label=f"trail {i} ({len(t)} edges)"))
    else:
        elements.append(
            Line2D([0], [0], color=TRAIL_COLORS[0], lw=2.0, label="Eulerian path/circuit")
        )
    if a.trails:
        for label in ("start", "end"):
            elements.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor=SURFACE,
                    markeredgecolor=INK,
                    markersize=9,
                    label=label,
                )
            )
    elements.append(Line2D([0], [0], color=FUTURE, lw=1.2, label="future edge"))

    fig.legend(
        handles=elements,
        loc="upper left",
        bbox_to_anchor=(0.02, 0.98),
        frameon=True,
        facecolor=SURFACE,
        edgecolor=RULE,
        fontsize=8.5,
    )
    fig.suptitle(
        f"{a.headline}\n{settings_text(a.settings)}",
        fontsize=12,
        fontweight="bold",
        color=INK,
        y=0.96,
        wrap=True,
    )
    note = " · ".join(a.notes)
    if note:
        fig.text(
            0.5,
            0.058,
            note.strip(" |"),
            ha="center",
            va="bottom",
            fontsize=8,
            color=MUTED,
        )


def draw_figure(
    fig: Figure,
    a: Analysis,
    layout: dict[str, tuple[float, float]] | None,
    upto: int | None = None,
) -> None:
    """Draw the walk graph with its Eulerian trail state on a matplotlib Figure."""
    fig.clear()
    ax = fig.add_subplot(111)
    ax.set_facecolor(SURFACE)
    ax.axis("off")
    if layout is None:
        ax.text(
            0.5,
            0.5,
            f"{a.headline}\n\n(The graph is not drawn: "
            f"{' and '.join(_counts_over_limit(_walk_graph(a)))}.\n"
            "Use the tables below to inspect the Eulerian trails and assembly contigs.)",
            ha="center",
            va="center",
            color=INK,
            fontsize=12,
            transform=ax.transAxes,
            gid="note",
        )
        return
    plot = _build_plot(ax, a, layout, FIGURE_ASPECT)
    plot.style(sum(len(trail) for trail in a.trails) if upto is None else upto)
    _legend_and_title(fig, a)


def save_figure(a: Analysis, path: str | os.PathLike) -> None:
    """Draw the graph and save it; the format comes from the extension."""
    extension = Path(path).suffix.lower().lstrip(".")
    if extension not in FIGURE_FORMATS:
        raise EuleraseError(
            f"Cannot save a figure to '{Path(path).name}'; use a file name ending in "
            + ", ".join(f".{f}" for f in FIGURE_FORMATS)
            + "."
        )
    fig = Figure(figsize=(10, 7))
    draw_figure(fig, a, compute_layout(a, FIGURE_ASPECT))
    fig.savefig(path, format=extension)

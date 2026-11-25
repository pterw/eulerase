"""Layout and figure of the walk graph: node coordinates and a matplotlib figure, no window."""

import math
import os
import textwrap
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.text import Text
from matplotlib.ticker import MaxNLocator

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


def trail_style(trail: int) -> tuple[str, str]:
    """The colour and the dash of trail `trail`, counted from 1: no two of the first nine trails
    share both."""
    colour = TRAIL_COLORS[(trail - 1) % len(TRAIL_COLORS)]
    return colour, TRAIL_DASHES[((trail - 1) // len(TRAIL_COLORS)) % len(TRAIL_DASHES)]


@dataclass(frozen=True)
class LegendEntry:
    """One line of the legend: its label, the colour and dash of its sample, and what it explains
    (`kind` is "trail", "start" or "end")."""

    label: str
    color: str
    dash: str
    kind: str


def legend_entries(a: Analysis) -> tuple[LegendEntry, ...]:
    """The legend of the figure and of the window's strip: each trail by its number, then the
    solid ring of a start and the dashed ring of an end."""
    trails = tuple(
        LegendEntry(f"trail {t}", *trail_style(t), "trail") for t in range(1, len(a.trails) + 1)
    )
    return trails + (
        LegendEntry("start ring", INK, TRAIL_DASHES[0], "start"),
        LegendEntry("end ring", INK, TRAIL_DASHES[1], "end"),
    )


def _counts_over_limit(a: Analysis) -> list[str]:
    """Each count of the drawing limit that `a` is over, with its number and its limit."""
    nodes = len(a.nodes)
    arrows = sum(len(t) for t in a.trails)
    over = []
    if nodes > MAX_DRAW_NODES:
        over.append(f"{nodes} nodes, over the limit of {MAX_DRAW_NODES}")
    if arrows > MAX_DRAW_EDGES:
        over.append(f"{arrows} arrows, over the limit of {MAX_DRAW_EDGES}")
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
    if _counts_over_limit(a):
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
BAND_WIDTH = 7.0  # line width of the lemon band under the current arrow
NODE_LINE = 1.5  # outline of a node
CURRENT_NODE_LINE = 3.0  # outline of the current node: wider, so the cue does not rely on colour
MARK_PITCH = 6.0  # the most room one copy mark takes along its line, in points
MARK_RADIUS = 0.4  # of a mark, in pitches
COUNT_PT = 5.5  # font size of the count of copies used, above the copy limit


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
    """An arrow with what is drawn for it: its lemon band, its label, and under the label either
    one mark for each copy or a count of the copies used. Everything exists from the start."""

    patch: FancyArrowPatch
    band: FancyArrowPatch
    label: Text | None
    marks: list[Circle]
    count: Text | None
    steps: tuple[tuple[int, int], ...]  # (global step, trail) of each copy
    weight: float

    def style(self, upto: int) -> None:
        """Restyle the arrow, its band, its label and its marks for the first `upto` steps."""
        used = [trail for step, trail in self.steps if step <= upto]
        active = any(step == upto for step, _ in self.steps)
        colour, line = trail_style(used[0]) if used else (FUTURE, TRAIL_DASHES[0])
        self.patch.set_color(colour)
        self.patch.set_linestyle(line)
        self.patch.set_linewidth(self.weight * (2.4 if active else 1.6 if used else 1.0))
        self.patch.set_mutation_scale(self.weight * (11 if used else 8))
        self.patch.set_zorder(3 if used else 1)
        self.band.set_visible(active)
        if self.label is not None:
            self.label.set_color(INK if used else MUTED)
            self.label.set_fontweight("bold" if active else "normal")
            self.label.get_bbox_patch().set_facecolor(LEMON if active else SURFACE)
        for copy, mark in enumerate(self.marks):
            mark.set_fill(copy < len(used))
        if self.count is not None:
            self.count.set_text(f"{len(used)}/{len(self.steps)}")


@dataclass
class _Plot:
    arrows: list[_Arrow]
    nodes: dict[str, Circle]
    targets: dict[int, str]  # the node each step of the walk arrives at, by global step
    overview: bool
    placement: Placement  # of the labels and words: what the marks go by, and the counts

    def style(self, upto: int) -> None:
        """Restyle the arrows, labels, marks and nodes for the first `upto` steps of the walk."""
        for arrow in self.arrows:
            arrow.style(upto)
        for node, disc in self.nodes.items():
            current = self.targets.get(upto) == node
            disc.set_facecolor(LEMON if current else INK if self.overview else SURFACE)
            disc.set_linewidth(CURRENT_NODE_LINE if current else 0 if self.overview else NODE_LINE)


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


def _add_node(ax, node: str, centre: Point, sizes: _Sizes, overview: bool) -> Circle:
    x, y = centre[0] / sizes.scale, centre[1] / sizes.scale
    disc = Circle(
        (x, y), sizes.radius / sizes.scale, edgecolor=INK, zorder=4, gid=f"node:{node}"
    )  # its fill and outline are the plot's `style`
    ax.add_patch(disc)
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
            parse_math=False,
        )
    return disc


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
        parse_math=False,
        bbox=None
        if is_word
        else dict(boxstyle="round,pad=0.15", facecolor=SURFACE, edgecolor="none", alpha=0.85),
    )


def _add_marks(ax, plan: _Plan, placed: Placed, scale: float) -> tuple[list[Circle], Text | None]:
    """The line under a label that kept room for it: one mark for each copy of the k-mer, or, above
    MAX_COPY_MARKS copies, one text for the copies used. It lies in the lowest MARKS_LINE_H of the
    label's box, inside its padding."""
    x0, y0, x1, _ = placed.box
    x, y = (x0 + x1) / 2 / scale, (y0 + LABEL_PAD + MARKS_LINE_H / 2) / scale
    copies = len(plan.steps)
    if copies > MAX_COPY_MARKS:
        widest = len(f"{copies}/{copies}") * MONO_ADVANCE
        size = min(COUNT_PT, (x1 - x0) / widest)  # the count stays inside the label's width
        count = ax.text(
            x, y, "", ha="center", va="center", family=MONO, fontsize=size, color=INK, zorder=5
        )
        count.set_gid(f"count:{plan.name}")
        return [], count
    pitch = min(MARK_PITCH, (x1 - x0) / copies)
    marks = [
        Circle(
            (x + (i - (copies - 1) / 2) * pitch / scale, y),
            MARK_RADIUS * pitch / scale,
            edgecolor=INK,
            facecolor=trail_style(plan.steps[0][1])[0],
            linewidth=0.8,
            zorder=5,
            gid=f"mark:{plan.name}",
        )
        for i in range(copies)
    ]
    for mark in marks:
        ax.add_patch(mark)
    return marks, None


def _add_arrow(
    ax, plan: _Plan, scale: float, label: Text | None, here: Placed | None, weight: float
) -> _Arrow:
    """The arrow of `plan` with its band under it and, when its label (`label`, placed at `here`)
    is drawn with a line of marks, those marks."""
    ends = ((plan.p0[0] / scale, plan.p0[1] / scale), (plan.p1[0] / scale, plan.p1[1] / scale))
    shape = dict(connectionstyle=f"arc3,rad={plan.rad}", shrinkA=0, shrinkB=0)
    patch = FancyArrowPatch(*ends, arrowstyle="-|>", gid=plan.gid, **shape)
    band = FancyArrowPatch(
        *ends,
        arrowstyle="-",
        color=LEMON,
        linewidth=BAND_WIDTH * weight,
        zorder=2,
        visible=False,
        gid=f"band:{plan.name}",
        **shape,
    )
    ax.add_patch(band)
    ax.add_patch(patch)
    marks, count = _add_marks(ax, plan, here, scale) if here and here.marks else ([], None)
    return _Arrow(patch, band, label, marks, count, plan.steps, weight)


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
        return _Plot([], {}, {}, overview, Placement({}, 0, 0))
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

    nodes = {node: _add_node(ax, node, centre, sizes, overview) for node, centre in pos.items()}
    for (kind, t, node), radius in zip(rings, ring_radii):
        _add_ring(ax, kind, t, pos[node], radius, scale)
    placed = placement.placed
    texts = {key: _add_text(ax, key, text, scale) for key, text in placed.items()}
    weight = OVERVIEW_WEIGHT if overview else 1.0
    arrows = [
        _add_arrow(ax, p, scale, texts.get(p.label_gid), placed.get(p.label_gid), weight)
        for p in plans
    ]
    targets = {s.global_step: s.target for trail in a.trails for s in trail}
    return _Plot(arrows, nodes, targets, overview, placement)


# ---------------------------------------------------------------------------
# The frame of a figure: headline, reading line and settings above the plot; the left-out text and
# the legend below it. Sizes come from character counts and font sizes, never from the renderer.
# ---------------------------------------------------------------------------
FRAME_MARGIN = 8.0  # points between the frame and the edge of the figure
FRAME_GAP = 3.0  # points between two parts of the frame
HEADLINE_PT, FRAME_PT = 10.5, 7.5  # font sizes of the headline and of every other frame text
LINE_PITCH = 1.3  # height of a line of frame text, in em, descenders and line gap included
SANS_ADVANCE, BOLD_ADVANCE = 0.58, 0.68  # widest likely average of a character, in em
LEGEND_HANDLE_EM = 2.8  # the sample and the room round it, in em
LEGEND_COLUMN_GAP_EM = 2.0
READING = "Grey arrows are not walked yet. One mark for each copy of a k-mer; filled when used."


def _wrapped(text: str, pt: float, advance: float, width: float) -> tuple[str, int]:
    """`text` broken into lines that fit `width` points by character count, and their number."""
    wrapped = textwrap.fill(text, width=max(8, int(width / (pt * advance))))
    return wrapped, wrapped.count("\n") + 1


def _draw_header(fig: Figure, a: Analysis, graph: bool) -> float:
    """The headline, the reading line (only over a drawn graph, `graph`) and the settings line,
    from the top left, wrapped to the figure's width. Returns the height they take, in points."""
    width, height = fig.get_figwidth() * 72, fig.get_figheight() * 72
    used = FRAME_MARGIN
    for gid, text, pt, bold in (
        ("headline", a.headline, HEADLINE_PT, True),
        ("reading", READING if graph else "", FRAME_PT, False),
        ("settings", settings_text(a.settings), FRAME_PT, False),
    ):
        if not text:
            continue
        advance = BOLD_ADVANCE if bold else SANS_ADVANCE
        lines, count = _wrapped(text, pt, advance, width - 2 * FRAME_MARGIN)
        fig.text(
            FRAME_MARGIN / width,
            1 - used / height,
            lines,
            ha="left",
            va="top",
            fontsize=pt,
            fontweight="bold" if bold else "normal",
            color=INK if bold else MUTED,
            gid=gid,
        )
        used += count * pt * LINE_PITCH + FRAME_GAP
    return used


def _legend_columns(labels: list[str], width: float) -> int:
    """How many columns of the legend fit `width` points, taking every column as wide as the
    widest label: at least 1, at most one for each entry."""
    unit = FRAME_PT * SANS_ADVANCE
    gap = LEGEND_COLUMN_GAP_EM * FRAME_PT
    sample = LEGEND_HANDLE_EM * FRAME_PT
    fit = int((width + gap) // (sample + gap + unit * max(map(len, labels))))
    return max(1, min(len(labels), fit))


def _draw_legend(fig: Figure, a: Analysis) -> float:
    """The legend of `legend_entries`, along the bottom edge. Returns the height it takes."""
    entries = legend_entries(a)
    width, height = fig.get_figwidth() * 72, fig.get_figheight() * 72
    columns = _legend_columns([e.label for e in entries], width - 2 * FRAME_MARGIN)
    rows = math.ceil(len(entries) / columns)
    handles = [
        Line2D([0], [0], color=e.color, linestyle=e.dash, lw=2.2, label=e.label) for e in entries
    ]
    fig.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(FRAME_MARGIN / width, FRAME_MARGIN / height),
        borderaxespad=0,
        ncols=columns,
        frameon=False,
        fontsize=FRAME_PT,
    )
    return rows * FRAME_PT * LINE_PITCH + (rows - 1) * 0.5 * FRAME_PT + 0.8 * FRAME_PT


def _draw_left_out(fig: Figure, plot: _Plot, bottom: float) -> None:
    """One text for the labels shortened and the labels left out, when either is above 0, with
    its lower edge `bottom` points above the figure's."""
    shortened, left_out = plot.placement.shortened, plot.placement.left_out
    if shortened or left_out:
        fig.text(
            FRAME_MARGIN / (fig.get_figwidth() * 72),
            bottom / (fig.get_figheight() * 72),
            f"Labels: {shortened} shortened (no copy marks), {left_out} left out.",
            ha="left",
            va="bottom",
            fontsize=FRAME_PT,
            color=MUTED,
            gid="left-out",
        )


def _draw_note(ax, a: Analysis, tables: bool) -> None:
    """The one explanation of the drawing limit, with both counts, in place of the plot. The window
    has tables and points at them (`tables`) and no header, so its note carries the headline; a
    saved figure has a header and no tables."""
    limits = " and ".join(_counts_over_limit(a))
    note = f"(The graph is not drawn: {limits}."
    if tables:
        hint = "Use the tables below to inspect the Eulerian trails and assembly contigs."
        note = f"{a.headline}\n\n{note}\n{hint}"
    ax.text(
        0.5,
        0.5,
        note + ")",
        ha="center",
        va="center",
        color=INK,
        fontsize=12,
        wrap=True,
        transform=ax.transAxes,
        gid="note",
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
    width, height = fig.get_figwidth() * 72, fig.get_figheight() * 72
    top = _draw_header(fig, a, graph=layout is not None)
    if layout is None:
        _draw_note(ax, a, tables=False)
        return
    footer = FRAME_MARGIN + _draw_legend(fig, a) + FRAME_GAP
    bottom = footer + FRAME_PT * LINE_PITCH + FRAME_GAP  # the left-out text lies above the legend
    ax.set_position(
        (
            FRAME_MARGIN / width,
            bottom / height,
            1 - 2 * FRAME_MARGIN / width,
            max(0.05, 1 - (top + bottom) / height),
        )
    )
    plot = _build_plot(ax, a, layout, FIGURE_ASPECT)
    plot.style(sum(len(trail) for trail in a.trails) if upto is None else upto)
    _draw_left_out(fig, plot, footer)


SPECTRUM_HEADROOM = 0.2  # of the tallest bar's height, added above it when a minimum is marked


def draw_spectrum(fig: Figure, a: Analysis) -> None:
    """Draw the k-mer spectrum: one bar per multiplicity, the minimum multiplicity marked.

    Bars below the minimum (the k-mers that setting removes) are grey and hatched; the rest are
    blue and plain, so colour is not the only cue.
    """
    fig.clear()
    fig.set_layout_engine("constrained")
    ax = fig.add_subplot(111)
    ax.set_facecolor(SURFACE)
    minimum = a.settings.min_multiplicity
    for times, count in a.spectrum:
        cut = minimum is not None and times < minimum
        ax.bar(
            times,
            count,
            color=FUTURE if cut else TRAIL_COLORS[0],
            edgecolor=INK,
            hatch="//" if cut else None,
        )
    if minimum is not None:
        ax.margins(y=SPECTRUM_HEADROOM)  # room above the tallest bar for the label
        line = minimum - 0.5
        ax.axvline(line, color=INK, linewidth=1.5)
        # the label lies along the top, on the side of the line that has more room, so it stays
        # inside the plot and clear of the y tick labels
        right = line < sum(ax.get_xlim()) / 2
        ax.text(
            line,
            0.98,  # a little below the top edge, so the label's box lies inside the axes
            f" minimum {minimum}" if right else f"minimum {minimum} ",
            color=INK,
            ha="left" if right else "right",
            va="top",
            transform=ax.get_xaxis_transform(),
        )
    ax.set_xlabel("Times seen", color=INK)
    ax.set_ylabel("Distinct k-mers", color=INK)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=1))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))


HOVER_PT = 8.0  # font size of the read-out of a node under the pointer


class StepScene:
    """The plot of `draw_figure`, drawn once on `fig` and restyled for each step: `show` changes
    properties of elements made here and creates none, moves no label and never lays out again.
    The headline, the legend and the other frame texts are not drawn: the window has its own."""

    def __init__(self, fig: Figure, a: Analysis, layout: dict[str, Point] | None) -> None:
        fig.clear()
        self._ax = ax = fig.add_axes((0, 0, 1, 1))
        ax.set_facecolor(SURFACE)
        ax.axis("off")
        self._plot: _Plot | None = None
        self._hover: Text | None = None
        if layout is None:
            _draw_note(ax, a, tables=True)
            return
        self._plot = _build_plot(ax, a, layout, fig.get_figwidth() / fig.get_figheight())
        self._plot.style(0)
        self._degrees = {row.node: (row.in_degree, row.out_degree) for row in a.nodes}
        self._hover = ax.text(
            0.01,
            0.99,
            "",
            transform=ax.transAxes,
            ha="left",
            va="top",
            family=MONO,
            fontsize=HOVER_PT,
            color=INK,
            zorder=7,
            visible=False,
            gid="hover",
            parse_math=False,
            bbox=dict(boxstyle="round,pad=0.3", facecolor=SURFACE, edgecolor=RULE),
        )

    def show(self, upto: int) -> None:
        """Style the plot for the first `upto` steps of the walk."""
        if self._plot is not None:
            self._plot.style(upto)

    @property
    def left_out(self) -> tuple[int, int]:
        """The labels shortened and the labels left out: the counts `draw_figure` writes."""
        if self._plot is None:
            return (0, 0)
        return (self._plot.placement.shortened, self._plot.placement.left_out)

    def node_at(self, x: float, y: float) -> str | None:
        """The node whose disc holds the display point (x, y), or None. Discs never overlap."""
        if self._plot is None:
            return None
        self._ax.apply_aspect()  # a figure resized and not yet drawn has a stale data transform
        to_display = self._ax.transData.transform
        for node, disc in self._plot.nodes.items():
            cx, cy = to_display(disc.center)
            reach = to_display((disc.center[0] + disc.radius, disc.center[1]))[0] - cx
            if math.dist((cx, cy), (x, y)) <= reach:
                return node
        return None

    def show_hover(self, node: str | None) -> None:
        """Show the read-out of `node` (its letters and its in and out degrees), or hide it."""
        if self._hover is None:
            return
        known = node in self._degrees
        if known:
            inward, outward = self._degrees[node]
            self._hover.set_text(f"{node}  in {inward}  out {outward}")
        self._hover.set_visible(known)


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

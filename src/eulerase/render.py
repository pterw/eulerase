"""Layout and figure of the walk graph: node coordinates and a matplotlib figure, no window."""

import math
import os
from collections import defaultdict
from pathlib import Path

import networkx as nx
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch

from eulerase.analysis import Analysis, settings_text
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


def _loop_geometry(x, y, r, i):
    tilt = 0.55 + 0.35 * i
    cx_loop = x + 2.2 * r * math.cos(tilt)
    cy_loop = y + 2.2 * r * math.sin(tilt)
    loop_r = 1.35 * r
    ang1 = tilt + math.pi - 0.7
    ang2 = tilt - math.pi + 0.7
    p_start = (cx_loop + loop_r * math.cos(ang1), cy_loop + loop_r * math.sin(ang1))
    p_end = (cx_loop + loop_r * math.cos(ang2), cy_loop + loop_r * math.sin(ang2))
    p_mid = (cx_loop + loop_r * math.cos(tilt), cy_loop + loop_r * math.sin(tilt))
    return p_start, p_end, p_mid, loop_r


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
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")

    g = _walk_graph(a)
    pos = layout
    if pos is None:
        ax.text(
            0.5,
            0.5,
            f"{a.headline}\n\n(The graph is not drawn: "
            f"{' and '.join(_counts_over_limit(g))}.\n"
            "Use the tables below to inspect the Eulerian trails and assembly contigs.)",
            ha="center",
            va="center",
            color=INK,
            fontsize=12,
            transform=ax.transAxes,
        )
        return

    step_of = {
        (s.source, s.target, s.key): (s.trail, s.step, s.global_step)
        for trail in a.trails
        for s in trail
    }
    starts = {t[0].source for t in a.trails}
    ends = {t[-1].target for t in a.trails}
    total = len(step_of)

    n_nodes = max(len(pos), 1)
    dists = [
        math.hypot(pos[v][0] - pos[u][0], pos[v][1] - pos[u][1])
        for u, v in g.edges()
        if u != v and u in pos and v in pos
    ]
    min_dist = min(dists) if dists else 0.25
    max_r_from_dist = min_dist * 0.38
    base_r = max(0.010, min(0.09, 0.55 / math.sqrt(n_nodes), max(0.010, max_r_from_dist)))
    font_s = max(5.5, min(9.5, 4.2 + base_r * 50))
    edge_fs = max(4.8, min(8.0, font_s - 1.0))
    multi_trail = len(a.trails) > 1
    upto_step = total if upto is None else upto

    cur_edge_key = None
    if upto_step > 0:
        for k_tuple, (_, _, g_step) in step_of.items():
            if g_step == upto_step:
                cur_edge_key = k_tuple
                break

    # Group parallel edges
    by_pair = defaultdict(list)
    for u, v, key in g.edges(keys=True):
        by_pair[(u, v)].append(key)

    rad_shifts = {}
    for (u, v), keys in by_pair.items():
        if u == v:
            continue
        has_rev = (v, u) in by_pair
        n_par = len(keys)
        base = 0.28 if has_rev else 0.0
        step = 0.24
        for i, k_id in enumerate(keys):
            rad_shifts[(u, v, k_id)] = (
                base + i * step if has_rev else -(n_par - 1) * step / 2 + i * step
            )

    # Draw edges
    for (u, v), keys in by_pair.items():
        if u == v:
            x, y = pos[u]
            for i, k_id in enumerate(keys):
                p_start, p_end, p_mid, loop_r = _loop_geometry(x, y, base_r, i)
                ti, _, g_step = step_of.get((u, v, k_id), (1, 1, 999999))
                past = g_step <= upto_step
                active = (u, v, k_id) == cur_edge_key
                color = (
                    (TRAIL_COLORS[(ti - 1) % len(TRAIL_COLORS)] if multi_trail else TRAIL_COLORS[0])
                    if past
                    else FUTURE
                )
                lw = 2.4 if active else 1.6 if past else 1.0
                patch = FancyArrowPatch(
                    p_start,
                    p_end,
                    connectionstyle="arc3,rad=-1.3",
                    arrowstyle="-|>",
                    mutation_scale=11 if past else 8,
                    color=color,
                    lw=lw,
                    zorder=3 if past else 1,
                )
                ax.add_patch(patch)
                kmer = g.edges[u, v, k_id]["kmer"]
                lbl = f"{kmer} (#{g_step})" if past else kmer
                ax.text(
                    p_mid[0],
                    p_mid[1] + 0.03,
                    lbl,
                    ha="center",
                    va="bottom",
                    fontsize=edge_fs,
                    color=INK if past else MUTED,
                    weight="bold" if active else "normal",
                    zorder=5,
                    bbox=dict(
                        boxstyle="round,pad=0.15", facecolor=SURFACE, edgecolor="none", alpha=0.85
                    ),
                )
            continue

        x1, y1 = pos[u]
        x2, y2 = pos[v]
        for k_id in keys:
            rad = rad_shifts.get((u, v, k_id), 0.0)
            ti, _, g_step = step_of.get((u, v, k_id), (1, 1, 999999))
            past = g_step <= upto_step
            active = (u, v, k_id) == cur_edge_key
            color = (
                (TRAIL_COLORS[(ti - 1) % len(TRAIL_COLORS)] if multi_trail else TRAIL_COLORS[0])
                if past
                else FUTURE
            )
            lw = 2.4 if active else 1.6 if past else 1.0
            patch = FancyArrowPatch(
                (x1, y1),
                (x2, y2),
                connectionstyle=f"arc3,rad={rad}",
                arrowstyle="-|>",
                mutation_scale=11 if past else 8,
                color=color,
                lw=lw,
                shrinkA=base_r * 80,
                shrinkB=base_r * 80,
                zorder=3 if past else 1,
            )
            ax.add_patch(patch)
            x12, y12 = (x1 + x2) / 2.0, (y1 + y2) / 2.0
            dx, dy = x2 - x1, y2 - y1
            dist = math.hypot(dx, dy) or 1.0
            mx = x12 + 0.5 * rad * dy
            my = y12 - 0.5 * rad * dx
            if abs(rad) < 1e-4:
                lx = mx - (dy / dist) * 0.04
                ly = my + (dx / dist) * 0.04
            else:
                lx = mx + (dy / dist) * (0.04 if rad > 0 else -0.04)
                ly = my - (dx / dist) * (0.04 if rad > 0 else -0.04)

            # In dense graphs (>25 edges), avoid plastering dozens or hundreds of
            # overlapping text boxes over straight segments or repeated copies.
            # Keep labels on the active step or self-loops in compact graphs.
            show_label = True
            if g.number_of_edges() > 25:
                show_label = active or (u == v and len(g) <= 40)

            if show_label:
                kmer = g.edges[u, v, k_id]["kmer"]
                lbl = f"{kmer} (#{g_step})" if past else kmer
                ax.text(
                    lx,
                    ly,
                    lbl,
                    ha="center",
                    va="center",
                    fontsize=edge_fs,
                    color=INK if past else MUTED,
                    weight="bold" if active else "normal",
                    zorder=5,
                    bbox=dict(
                        boxstyle="round,pad=0.15", facecolor=SURFACE, edgecolor="none", alpha=0.85
                    ),
                )

    # Identify landmark branch / junction nodes
    branch_nodes = set()
    for n in g:
        if len(set(g.predecessors(n))) > 1 or len(set(g.successors(n))) > 1:
            branch_nodes.add(n)

    # Draw nodes
    for node, (x, y) in pos.items():
        is_start = node in starts
        is_end = node in ends
        is_branch = node in branch_nodes
        circ = Circle(
            (x, y),
            base_r * (1.3 if (is_start or is_end or is_branch) and n_nodes > 40 else 1.0),
            facecolor=SURFACE,
            edgecolor=INK,
            linewidth=1.2 if (is_start or is_end or is_branch) and n_nodes > 40 else 1.5,
            zorder=4,
        )
        ax.add_patch(circ)

        # In dense graphs (>40 nodes), suppress intermediate non-branching chain labels
        # to avoid solid walls of text; preserve start, end, and multi-way branch nodes.
        should_label_node = True
        if n_nodes > 40:
            should_label_node = is_start or is_end or is_branch

        if should_label_node:
            lbl = (
                node if n_nodes <= 40 else (node[:3] + ".." + node[-2:] if len(node) > 7 else node)
            )
            ax.text(
                x,
                y,
                lbl,
                ha="center",
                va="center",
                fontsize=font_s,
                weight="bold",
                color=INK,
                zorder=6,
                bbox=dict(
                    boxstyle="round,pad=0.15", facecolor=SURFACE, edgecolor=INK, lw=0.6, alpha=0.9
                )
                if n_nodes > 40
                else None,
            )

    ax.margins(0.12)
    ax.autoscale_view()

    # Legend / Title
    elements = []
    if multi_trail:
        for i, t in enumerate(a.trails, 1):
            c = TRAIL_COLORS[(i - 1) % len(TRAIL_COLORS)]
            elements.append(Line2D([0], [0], color=c, lw=2.2, label=f"trail {i} ({len(t)} edges)"))
    else:
        elements.append(
            Line2D([0], [0], color=TRAIL_COLORS[0], lw=2.0, label="Eulerian path/circuit")
        )
    if any(node in starts for node in pos):
        elements.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=SURFACE,
                markeredgecolor=INK,
                markersize=9,
                label="start",
            )
        )
    if any(node in ends for node in pos):
        elements.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=SURFACE,
                markeredgecolor=INK,
                markersize=9,
                label="end",
            )
        )
    if any(node in branch_nodes for node in pos) and n_nodes > 40:
        elements.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=SURFACE,
                markeredgecolor=INK,
                markersize=9,
                label="branch / junction",
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

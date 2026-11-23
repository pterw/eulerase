"""Layout and figure of the walk graph: node coordinates and a matplotlib figure, no window."""

import math
import os
import random
from collections import defaultdict
from pathlib import Path

import networkx as nx
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch

from eulerase.analysis import Analysis, settings_text
from eulerase.errors import EuleraseError

# ---------------------------------------------------------------------------
# Palette (citrus) and limits
# ---------------------------------------------------------------------------
PITH = "#FFF8E1"  # background
LEMON = "#F9D923"  # ordinary node
RIND = "#E3B505"  # node outline
LIME = "#6BA32B"  # trail start
ORANGE = "#F28C28"  # trail end
EDGE = "#8C6D1F"  # edges (single trail)
INK = "#3B2F0B"  # text
FUTURE = "#D6CBA0"  # not-yet-traversed edges
FUTURE_TXT = "#B3A678"
PANEL = "#FFF1BF"
TRAIL_COLORS = ["#2F7D32", "#D9480F", "#C2255C", "#1C7ED6", "#9C36B5", "#8C6D1F"]
BRANCH = "#B53A55"  # branch / junction node in dense graphs

MAX_DRAW_NODES = 400  # above this the graph is not drawn (tables still work)
FIGURE_FORMATS = ("png", "svg", "pdf")


def _walk_graph(a: Analysis) -> nx.MultiDiGraph:
    """The walk graph again, from the nodes and steps the analysis holds."""
    g = nx.MultiDiGraph()
    g.add_nodes_from(row.node for row in a.nodes)
    for trail in a.trails:
        for s in trail:
            g.add_edge(s.source, s.target, key=s.key, kmer=s.kmer)
    return g


def _rotate_horizontal(pos: dict[str, tuple[float, float]]):
    """Align the principal axis of the layout with the horizontal axis."""
    pts = list(pos.values())
    if len(pts) < 2:
        return pos
    mean_x = sum(x for x, y in pts) / len(pts)
    mean_y = sum(y for x, y in pts) / len(pts)
    sxx = sum((x - mean_x) ** 2 for x, y in pts)
    syy = sum((y - mean_y) ** 2 for x, y in pts)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in pts)
    angle = 0.5 * math.atan2(2 * sxy, sxx - syy)
    cos_a, sin_a = math.cos(-angle), math.sin(-angle)
    return {
        k: (
            (x - mean_x) * cos_a - (y - mean_y) * sin_a,
            (x - mean_x) * sin_a + (y - mean_y) * cos_a,
        )
        for k, (x, y) in pos.items()
    }


def _component_layout(sub: nx.MultiDiGraph) -> dict[str, tuple[float, float]]:
    n = len(sub)
    if n == 1:
        return {next(iter(sub)): (0.0, 0.0)}
    if n == 2:
        nodes = list(sub)
        return {nodes[0]: (-1.0, 0.0), nodes[1]: (1.0, 0.0)}
    simple = nx.Graph(sub)
    seed = random.Random(2026).randint(0, 10000)
    pos = nx.spring_layout(
        simple,
        k=1.8 / math.sqrt(n),
        iterations=160,
        seed=seed,
        scale=2.0,
    )
    pos = _rotate_horizontal(pos)
    pts = list(pos.values())
    min_x, max_x = min(x for x, y in pts), max(x for x, y in pts)
    min_y, max_y = min(y for x, y in pts), max(y for x, y in pts)
    rx, ry = max(max_x - min_x, 1e-6), max(max_y - min_y, 1e-6)
    return {k: ((x - min_x) / rx, (y - min_y) / ry) for k, (x, y) in pos.items()}


def _ordered_subgraph(g: nx.MultiDiGraph, nodes: list[str]) -> nx.MultiDiGraph:
    """The subgraph on `nodes` (already sorted), built so its node and edge order is fixed."""
    sub = nx.MultiDiGraph()
    sub.add_nodes_from(nodes)
    for u in nodes:
        sub.add_edges_from(g.out_edges(u, keys=True, data=True))
    return sub


def compute_layout(a: Analysis) -> dict[str, tuple[float, float]] | None:
    """Layout components in horizontal rows."""
    g = _walk_graph(a)
    if len(g) > MAX_DRAW_NODES:
        return None
    # A component is a set, so its iteration order follows string hashing: take every node and
    # edge from the walk graph in sorted-node order, and order equal-size components by name.
    ordered = sorted(
        (sorted(c) for c in nx.weakly_connected_components(g)),
        key=lambda nodes: (-len(nodes), nodes[0]),
    )
    comps = [_ordered_subgraph(g, nodes) for nodes in ordered]
    if not comps:
        return {}
    layouts = [_component_layout(c) for c in comps]

    row_w = 0.0
    gap = 0.4
    placements = []
    for lyt in layouts:
        pts = list(lyt.values())
        w = max(x for x, y in pts) - min(x for x, y in pts) if len(pts) > 1 else 0.5
        placements.append((row_w, 0.0, lyt))
        row_w += max(w, 0.5) + gap

    out = {}
    for ox, oy, lyt in placements:
        for k, (x, y) in lyt.items():
            out[k] = (ox + x, oy + y)

    xs, ys = [p[0] for p in out.values()], [p[1] for p in out.values()]
    cx, cy = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2
    ext = max(max(xs) - min(xs), max(ys) - min(ys), 1e-9)
    f = 2.0 / ext
    return {k: ((px - cx) * f, (py - cy) * f) for k, (px, py) in out.items()}


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
    ax.set_facecolor(PITH)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")

    g = _walk_graph(a)
    pos = layout
    if pos is None:
        ax.text(
            0.5,
            0.5,
            f"{a.headline}\n\n(Graph has {len(g)} nodes; "
            f"drawing is disabled above {MAX_DRAW_NODES}.\n"
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
                    (TRAIL_COLORS[(ti - 1) % len(TRAIL_COLORS)] if multi_trail else EDGE)
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
                    color=INK if past else FUTURE_TXT,
                    weight="bold" if active else "normal",
                    zorder=5,
                    bbox=dict(
                        boxstyle="round,pad=0.15", facecolor=PITH, edgecolor="none", alpha=0.85
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
                (TRAIL_COLORS[(ti - 1) % len(TRAIL_COLORS)] if multi_trail else EDGE)
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
                    color=INK if past else FUTURE_TXT,
                    weight="bold" if active else "normal",
                    zorder=5,
                    bbox=dict(
                        boxstyle="round,pad=0.15", facecolor=PITH, edgecolor="none", alpha=0.85
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
        fill = (
            LIME
            if (is_start and not is_end)
            else ORANGE
            if (is_end and not is_start)
            else (BRANCH if is_branch and n_nodes > 40 else LEMON)
        )
        circ = Circle(
            (x, y),
            base_r * (1.3 if (is_start or is_end or is_branch) and n_nodes > 40 else 1.0),
            facecolor=fill,
            edgecolor=RIND,
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
                    boxstyle="round,pad=0.15", facecolor=PITH, edgecolor=RIND, lw=0.6, alpha=0.9
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
        elements.append(Line2D([0], [0], color=EDGE, lw=2.0, label="Eulerian path/circuit"))
    if any(node in starts for node in pos):
        elements.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=LIME,
                markeredgecolor=RIND,
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
                markerfacecolor=ORANGE,
                markeredgecolor=RIND,
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
                markerfacecolor=BRANCH,
                markeredgecolor=RIND,
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
        facecolor=PITH,
        edgecolor=PANEL,
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
            color=FUTURE_TXT,
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
    draw_figure(fig, a, compute_layout(a))
    fig.savefig(path, format=extension)

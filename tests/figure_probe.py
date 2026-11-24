"""Helpers the figure tests share: an Agg figure, elements found by `gid` role, extents.

Flat import, as `genomes.py`: `from figure_probe import ...`.
"""

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from eulerase.analysis import Settings, analyze
from eulerase.render import FIGURE_ASPECT, compute_layout, draw_figure
from eulerase.sequence import Record

SAVED_SIZE = (10, 7)  # inches: the size of a saved figure
SMALL_SIZE = (5.5, 3.4)  # inches: the plot area keeps its aspect, so it shrinks in pixels


def _reads(*reads: str, k: int):
    return analyze([Record(f"r{i}", text) for i, text in enumerate(reads, 1)], Settings(k=k))


# The fixed inputs: each name maps to a function so that a test builds only what it uses.
INPUTS = {
    "lemon": lambda: _reads("LEMONLEMONADE", k=3),
    "repeat": lambda: _reads("ACGTTGCAACGTTGCA", k=4),
    "tandem": lambda: _reads("AAAAAC", k=3),
    "overlap": lambda: _reads("ACGTGCAATG", "CAATGCCGTA", "CCGTAACGTG", k=5),
    # three trails of 12, 5 and 8 steps, none of which ends on the node it starts on
    "three": lambda: _reads("ACGTACGGA", "TTGCCATTAGGC", "GGATCCTTAAGCATGC", k=5),
    # one trail that ends where it starts
    "circuit": lambda: _reads("ACGTAC", k=3),
}


def drawn(a, size=SAVED_SIZE, upto=None, layout=None):
    """A figure drawn by `draw_figure` on an Agg canvas, rendered once. Returns (fig, layout)."""
    if layout is None:
        layout = compute_layout(a, FIGURE_ASPECT)
    fig = Figure(figsize=size)
    FigureCanvasAgg(fig)
    draw_figure(fig, a, layout, upto)
    fig.canvas.draw()
    return fig, layout


def elements(fig, role):
    """Every artist of the plot axes whose `gid` is `role` or starts with `role:`."""
    return [
        artist
        for artist in fig.axes[0].get_children()
        if (gid := artist.get_gid()) is not None and (gid == role or gid.startswith(role + ":"))
    ]


def extent(fig, artist):
    """The artist's box in display pixels, as (x0, y0, x1, y1)."""
    box = artist.get_window_extent(fig.canvas.get_renderer())
    return (box.x0, box.y0, box.x1, box.y1)


def overlap(first, second, margin=1.0):
    """True when two (x0, y0, x1, y1) boxes intersect or lie within `margin` pixels."""
    return (
        first[0] - margin < second[2]
        and second[0] - margin < first[2]
        and first[1] - margin < second[3]
        and second[1] - margin < first[3]
    )


def inside(box, fig):
    """True when the box lies inside the figure's own box."""
    return box[0] >= 0 and box[1] >= 0 and box[2] <= fig.bbox.width and box[3] <= fig.bbox.height

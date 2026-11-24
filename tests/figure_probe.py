"""Helpers the figure tests share: an Agg figure, elements found by `gid` role, extents.

Flat import, as `genomes.py`: `from figure_probe import ...`.
"""

from genomes import random_dna
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.text import Text

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
    # one trail in which the 3-mer ABC is seen three times
    "repeated": lambda: _reads("ABCXABCYABCZ", k=3),
    # three trails of 2, 11 and 6 steps; the 3-mer BAB is seen more than once in two of them
    "shared": lambda: _reads("AACBABABC", "BABABXC", "XABBBCBCB", k=3),
}


# Ten disjoint trails of ten different lengths (111 nodes: above the overview limit, which a test
# of the full drawing raises), so that the trail colour and the dash pattern both wrap.
def ten():
    reads = (random_dna(8 + 2 * i if i < 5 else 9 + 2 * i, seed=110 + i) for i in range(10))
    return _reads(*reads, k=5)


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


DASH_NAMES = {"-": "solid", "--": "dashed", ":": "dotted"}
DASH_NAMES |= {name: name for name in ("solid", "dashed", "dotted")}


def dash(artist):
    """The line style of a patch or a line, as one of "solid", "dashed", "dotted"."""
    return DASH_NAMES[artist.get_linestyle()]


def built(a, size=SAVED_SIZE, upto=None, layout=None):
    """Like `drawn`, without rendering: enough for the properties of the elements."""
    if layout is None:
        layout = compute_layout(a, FIGURE_ASPECT)
    fig = Figure(figsize=size)
    FigureCanvasAgg(fig)
    draw_figure(fig, a, layout, upto)
    return fig


def frame(fig, role):
    """The figure-level texts (outside the plot axes) whose `gid` is `role`."""
    return [t for t in fig.texts if t.get_gid() == role]


def legends(fig):
    """Every visible legend of the figure and of its axes."""
    found = fig.findobj(Legend)
    return [legend for legend in found if legend.get_visible()]


def visible_texts(fig):
    """Every visible, non-empty text of the figure, but not the legend's own label texts."""
    own = {text for legend in legends(fig) for text in legend.get_texts()}
    return [
        t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip() and t not in own
    ]

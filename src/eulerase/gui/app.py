"""The window frame: the form, the status bar and the busy state. It decides nothing.

Nothing here runs at import: no widget, no root, no thread. A `Window` is built from a root and
reaches its presenter only through `bind`, so no callback is wired before then.
"""

import tkinter
import tkinter.filedialog
import tkinter.font
import tkinter.messagebox
import tkinter.ttk as ttk
from dataclasses import dataclass
from pathlib import Path

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from eulerase import render
from eulerase.gui.presenter import (
    EXAMPLES,
    EXPORT_KINDS,
    MIN_WINDOW,
    FormValues,
    ResultView,
    TableView,
    fit_window,
)

MODES = (("sequences", "Sequences or FASTA"), ("kmers", "K-mer list"))
MIN_MULTIPLICITY_HELP = (
    "Drops k-mers seen fewer times before the contigs are built. For a set of reads this "
    "removes sequencing errors; on one assembled sequence a low count marks a rarer k-mer, "
    "not an error."
)
CLIP_TIPS_HELP = (
    "Removes a short dead-end branch (under 2k k-mers) when a better-supported branch leaves or "
    "enters the same node, before the contigs are built. For reads this can remove errors; it can "
    "also remove the true end of a sequence. It removes dead ends only."
)
SIDEBAR_TEXT_WIDTH = 250  # pixels at which the help sentences wrap
VIEW_BUTTONS = ("Zoom", "Pan", "Reset view")
STEP_BUTTONS = ("First", "Back", "Play", "Forward", "Last")  # the third reads "Pause" while playing
STEP_METHODS = ("step_start", "step_back", "toggle_play", "step_forward", "step_end")
ASSEMBLY = "Assembly"
WIDTH_ROWS = 50  # rows of a table read to size its columns
MAX_COLUMN_PX = 320  # a longer cell scrolls instead of widening the table
CELL_PAD_PX = 16
TK_KEY_NAMES = {"Ctrl": "Control", "Esc": "Escape", "Enter": "Return"}
DASH_PATTERNS = {"solid": (), "dashed": (6, 3), "dotted": (2, 3)}  # Tk dash, by line style name
MIN_CANVAS_PX = 2  # a canvas not yet mapped reports 1
RESIZE_DELAY_MS = 150  # a drag of the window edge lays out once, after it stops
ZOOM_IN, ZOOM_OUT = 0.85, 1.18  # one wheel notch, as a factor of the visible range


@dataclass(frozen=True)
class Command:
    """A menu entry: `method` names a Presenter method, called with `args`. `shortcut` is the
    text of the accelerator column ("Ctrl+Shift+F"), or None."""

    menu: str
    label: str
    shortcut: str | None
    method: str
    args: tuple = ()


def _export(kind: str, label: str, letter: str) -> Command:
    return Command("File", label, f"Ctrl+Shift+{letter}", "export", (kind,))


COMMANDS = (
    Command("File", "Open file...", "Ctrl+Shift+O", "open_file"),
    _export("figure", "Save figure...", "F"),
    _export("csv", "Save tables as CSV...", "C"),
    _export("fasta", "Save contigs as FASTA...", "A"),
    _export("genbank", "Save contigs as GenBank...", "B"),
    _export("gfa", "Save graph as GFA...", "G"),
    Command("File", "Quit", "Ctrl+Q", "quit"),
    *(Command("Examples", e.name, None, "load_example", (e.name,)) for e in EXAMPLES),
    Command("Walk", "Build", "Ctrl+Enter", "build"),
    Command("Walk", "Stop", "Esc", "stop"),
    Command("Walk", "Play or pause", "F5", "toggle_play"),
    Command("Walk", "Step back", "Alt+Left", "step_back"),
    Command("Walk", "Step forward", "Alt+Right", "step_forward"),
    Command("Walk", "First step", "Alt+Home", "step_start"),
    Command("Walk", "Last step", "Alt+End", "step_end"),
)


def tk_sequence(shortcut: str) -> str:
    """The Tk event sequence of a shortcut's text. Tk keysyms are case sensitive: a letter
    arrives in lower case unless Shift is held."""
    parts = shortcut.split("+")
    keys = [TK_KEY_NAMES.get(p, p) for p in parts]
    if "Shift" not in parts and len(keys[-1]) == 1:
        keys[-1] = keys[-1].lower()
    return "<" + "-".join(keys) + ">"


def _wrap_to_width(event) -> None:
    event.widget.configure(wraplength=max(event.width, 100))


def _column_width(font, header: str, cells) -> int:
    """The widest of the header and the cells given, plus padding, capped."""
    return min(max(font.measure(text) for text in (header, *cells)) + CELL_PAD_PX, MAX_COLUMN_PX)


class Window:
    """The frame of the window; `bind` attaches every callback to the presenter."""

    def __init__(self, root):
        self.root = root
        self.presenter = None
        root.title("Eulerase")
        width, height = fit_window(root.winfo_screenwidth(), root.winfo_screenheight())
        root.geometry(f"{width}x{height}")
        root.minsize(*MIN_WINDOW)
        root.configure(background=render.WINDOW)
        root.columnconfigure(1, weight=1)
        root.rowconfigure(0, weight=1)
        self._style()
        self._sidebar()
        self._explorer()
        self._graph()
        self._step_controls()
        self._status_bar()

    def _style(self) -> None:
        self.bold = tkinter.font.nametofont("TkDefaultFont").copy()
        self.bold.configure(weight="bold")
        heading_font = tkinter.font.nametofont("TkHeadingFont")
        style = ttk.Style(self.root)
        style.theme_use("clam")
        style.configure(
            ".", background=render.WINDOW, foreground=render.INK, bordercolor=render.RULE
        )
        style.configure("TLabel", background=render.WINDOW, foreground=render.INK)
        style.configure("Muted.TLabel", foreground=render.MUTED)
        style.configure("Error.TLabel", foreground=render.ERROR)
        style.configure("Heading.TLabel", font=self.bold)
        style.configure("TCheckbutton", background=render.WINDOW)
        style.configure("TRadiobutton", background=render.WINDOW)
        style.configure("TButton", padding=5, bordercolor=render.FIELD_BORDER)
        style.configure("TEntry", fieldbackground=render.SURFACE, bordercolor=render.FIELD_BORDER)
        style.configure(
            "Primary.TButton",
            background=render.LEMON,
            bordercolor=render.PRIMARY_BORDER,
            font=self.bold,
        )
        style.map("Primary.TButton", background=[("active", render.LEMON)])
        style.configure("Treeview", background=render.SURFACE, fieldbackground=render.SURFACE)
        style.configure("Treeview.Heading", background=render.WINDOW, font=heading_font)
        style.map(
            "Treeview",
            background=[("selected", render.LEMON_SOFT)],
            foreground=[("selected", render.INK)],
        )
        style.configure("TProgressbar", background=render.MUTED, troughcolor=render.RULE)

    def _sidebar(self) -> None:
        side = ttk.Frame(self.root, padding=12)
        side.grid(row=0, column=0, sticky="ns")
        self.mode = tkinter.StringVar(value=MODES[0][0])
        self.k_text = tkinter.StringVar(value="5")
        self.both_strands = tkinter.BooleanVar(value=False)
        self.unique = tkinter.BooleanVar(value=False)
        self.min_multiplicity_text = tkinter.StringVar(value="")
        self.clip_tips = tkinter.BooleanVar(value=False)

        ttk.Label(side, text="Input", style="Heading.TLabel").pack(anchor="w")
        self.text = tkinter.Text(
            side,
            width=34,
            height=11,
            wrap="char",
            font="TkFixedFont",
            background=render.SURFACE,
            foreground=render.INK,
            insertbackground=render.INK,
            relief="flat",
            highlightthickness=1,
            highlightbackground=render.FIELD_BORDER,
            highlightcolor=render.INK,
            undo=True,
        )
        self.text.pack(fill="x", pady=(4, 6))
        self.mode_buttons = [
            ttk.Radiobutton(side, text=label, variable=self.mode, value=value)
            for value, label in MODES
        ]
        for button in self.mode_buttons:
            button.pack(anchor="w")
        files = ttk.Frame(side)
        files.pack(fill="x", pady=8)
        self.open_button = ttk.Button(files, text="Open file")
        self.examples_button = ttk.Button(files, text="Examples")
        self.open_button.pack(side="left")
        self.examples_button.pack(side="left", padx=(8, 0))

        ttk.Label(side, text="Settings", style="Heading.TLabel").pack(anchor="w", pady=(8, 0))
        k_row = ttk.Frame(side)
        k_row.pack(fill="x", pady=4)
        ttk.Label(k_row, text="k (k-mer length)").pack(side="left")
        self.k_entry = ttk.Entry(k_row, textvariable=self.k_text, width=5)
        self.k_entry.pack(side="right")
        ttk.Checkbutton(side, text="Both strands (DNA only)", variable=self.both_strands).pack(
            anchor="w", pady=2
        )
        ttk.Checkbutton(side, text="Treat as a set (each k-mer once)", variable=self.unique).pack(
            anchor="w", pady=2
        )
        min_row = ttk.Frame(side)
        min_row.pack(fill="x", pady=(6, 0))
        ttk.Label(min_row, text="Minimum multiplicity").pack(side="left")
        ttk.Entry(min_row, textvariable=self.min_multiplicity_text, width=5).pack(side="right")
        self._help(side, MIN_MULTIPLICITY_HELP)
        ttk.Checkbutton(side, text="Clip tips", variable=self.clip_tips).pack(
            anchor="w", pady=(6, 0)
        )
        self._help(side, CLIP_TIPS_HELP)

        buttons = ttk.Frame(side)
        buttons.pack(side="bottom", fill="x", pady=(12, 0))
        self.build_button = ttk.Button(buttons, text="Build  Ctrl+Enter", style="Primary.TButton")
        self.stop_button = ttk.Button(buttons, text="Stop", state="disabled")
        self.build_button.pack(side="left")
        self.stop_button.pack(side="left", padx=(8, 0))

    @staticmethod
    def _help(parent, text: str) -> None:
        ttk.Label(parent, text=text, style="Muted.TLabel", wraplength=SIDEBAR_TEXT_WIDTH).pack(
            anchor="w", pady=(2, 0)
        )

    def _explorer(self) -> None:
        """The graph, its strip, the step controls and the tables: one frame, so a comparison
        view can later become a sibling tab. Later tasks fill the empty parts."""
        self.explorer = ttk.Frame(self.root, name="explorer", padding=(0, 8, 8, 0))
        self.explorer.grid(row=0, column=1, sticky="nsew")
        self.headline = ttk.Label(self.explorer, text="", style="Heading.TLabel")
        self.headline.pack(fill="x", padx=8, pady=(0, 6))
        divider = ttk.PanedWindow(self.explorer, orient="vertical")
        divider.pack(fill="both", expand=True)
        graph = ttk.Frame(divider)
        self.canvas_frame = tkinter.Frame(graph, background=render.SURFACE)
        self.canvas_frame.pack(fill="both", expand=True)
        self.strip = ttk.Frame(graph)
        self.strip.pack(fill="x", pady=4)
        self.legend_frame = ttk.Frame(self.strip)
        self.legend_frame.pack(side="left")
        self.view_buttons = {
            name: ttk.Button(self.strip, text=name, state="disabled") for name in VIEW_BUTTONS
        }
        for button in reversed(list(self.view_buttons.values())):
            button.pack(side="right", padx=(6, 0))
        self.steps_frame = ttk.Frame(graph)
        self.steps_frame.pack(fill="x")
        self.tables_frame = ttk.Frame(divider)
        divider.add(graph, weight=3)
        divider.add(self.tables_frame, weight=1)
        self.notebook = None  # made by the first result: its tabs are the tables of that result
        self.trees = {}
        self.notes = {}

    def _status_bar(self) -> None:
        bar = ttk.Frame(self.root, padding=(8, 2))
        bar.grid(row=1, column=0, columnspan=2, sticky="ew")
        self.progress = ttk.Progressbar(bar, mode="indeterminate", length=80)
        self.status = ttk.Label(bar, text="", justify="left")
        self.status.pack(side="left", fill="x", expand=True)
        bar.bind("<Configure>", self._wrap_status)

    def _wrap_status(self, event) -> None:
        self.status.configure(wraplength=max(event.width - 120, 100))

    def bind(self, presenter) -> None:
        """Store the presenter, then attach every callback of the frame to it."""
        self.presenter = presenter
        self.build_button.configure(command=presenter.build)
        self.stop_button.configure(command=presenter.stop)
        for button in self.mode_buttons:
            button.configure(command=lambda: presenter.mode_changed(self.mode.get()))
        self.text.bind("<Control-Return>", self._ctrl_enter)
        self.root.protocol("WM_DELETE_WINDOW", presenter.quit)
        self._menu_bar()
        open_command = next(c for c in COMMANDS if c.method == "open_file")
        self.open_button.configure(command=lambda: self._run(open_command))
        self.examples_button.configure(command=self._post_examples)
        for c in COMMANDS:
            if c.shortcut and c.method != "build":  # Ctrl+Enter is bound on the input box
                self.root.bind(tk_sequence(c.shortcut), lambda event, c=c: self._key(c))
        self._bind_graph(presenter)

    def _bind_graph(self, presenter) -> None:
        """The step controls call the presenter; the view controls and canvas are the window's."""
        for button, method in zip(self.step_buttons, STEP_METHODS, strict=True):
            button.configure(command=getattr(presenter, method))
        self.scale.configure(command=lambda value: presenter.set_step(round(float(value))))
        self.view_buttons["Zoom"].configure(command=self.toolbar.zoom)
        self.view_buttons["Pan"].configure(command=self.toolbar.pan)
        self.view_buttons["Reset view"].configure(command=self._reset_view)
        self.canvas_widget.bind("<Configure>", self._on_configure, add="+")  # keep matplotlib's
        self.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.canvas.mpl_connect("motion_notify_event", self._on_motion)

    def _menu_bar(self) -> None:
        bar = tkinter.Menu(self.root)
        menus = {}
        self._export_entries = {}
        for c in COMMANDS:
            if c.menu not in menus:
                menus[c.menu] = tkinter.Menu(bar, tearoff=False)
                bar.add_cascade(label=c.menu, menu=menus[c.menu])
            menu = menus[c.menu]
            menu.add_command(
                label=c.label, accelerator=c.shortcut or "", command=lambda c=c: self._run(c)
            )
            if c.method == "export":
                self._export_entries[c.args[0]] = (menu, menu.index("end"))
        self.examples_menu = menus["Examples"]
        self.root.configure(menu=bar)
        self.set_exports_enabled(dict.fromkeys(EXPORT_KINDS, False))

    def _run(self, command: Command) -> None:
        getattr(self.presenter, command.method)(*command.args)

    def _key(self, command: Command) -> str:
        self._run(command)
        return "break"

    def _post_examples(self) -> None:
        button = self.examples_button
        self.examples_menu.tk_popup(
            button.winfo_rootx(), button.winfo_rooty() + button.winfo_height()
        )

    def _ctrl_enter(self, event) -> str:
        self.presenter.build()
        return "break"  # no newline in the input box

    def close(self) -> None:
        self.root.destroy()

    def read_form(self) -> FormValues:
        return FormValues(
            text=self.text.get("1.0", "end-1c"),
            input_mode=self.mode.get(),
            k_text=self.k_text.get(),
            both_strands=self.both_strands.get(),
            unique=self.unique.get(),
            min_multiplicity_text=self.min_multiplicity_text.get(),
            clip_tips=self.clip_tips.get(),
        )

    def set_form(self, form: FormValues) -> None:
        self.text.delete("1.0", "end")
        self.text.insert("1.0", form.text)
        self.mode.set(form.input_mode)
        self.k_text.set(form.k_text)
        self.both_strands.set(form.both_strands)
        self.unique.set(form.unique)
        self.min_multiplicity_text.set(form.min_multiplicity_text)
        self.clip_tips.set(form.clip_tips)

    def set_k_enabled(self, enabled: bool) -> None:
        self.k_entry.state(["!disabled" if enabled else "disabled"])

    def show_status(self, message: str, error: bool = False) -> None:
        self.status.configure(text=message, style="Error.TLabel" if error else "TLabel")

    def set_exports_enabled(self, enabled) -> None:
        for kind, (menu, index) in self._export_entries.items():
            menu.entryconfigure(index, state="normal" if enabled[kind] else "disabled")

    def confirm(self, message: str) -> bool:
        return tkinter.messagebox.askyesno("Eulerase", message, parent=self.root)

    def ask_open_path(self) -> str | None:
        return tkinter.filedialog.askopenfilename(parent=self.root) or None

    def ask_save_path(self, kind: str, suggested: str) -> str | None:
        # The CSV name is a prefix: no file of that name is written, so no overwrite prompt.
        return (
            tkinter.filedialog.asksaveasfilename(
                parent=self.root,
                initialfile=suggested,
                defaultextension=Path(suggested).suffix,
                confirmoverwrite=kind != "csv",
            )
            or None
        )

    def confirm_export(self, paths, existing) -> bool:
        message = "Write these files?\n" + "\n".join(paths)
        if existing:
            message += "\n\nThese already exist and will be replaced:\n" + "\n".join(existing)
        return tkinter.messagebox.askyesno("Eulerase", message, parent=self.root)

    def save_figure(self, analysis, path: str) -> None:
        render.save_figure(analysis, path)  # EuleraseError or OSError: the presenter shows both

    def show_result(self, result: ResultView) -> None:
        self._show_graph(result)
        if self.notebook is None:
            self._make_tabs(result.tables)
        for table in result.tables:
            self._fill_table(table)
        self._set_details(result.details)
        render.draw_spectrum(self.spectrum, result.analysis)
        self.spectrum_canvas.draw_idle()

    def clear_result(self) -> None:
        self._clear_graph()
        if self.notebook is None:
            return
        for name, tree in self.trees.items():
            tree.delete(*tree.get_children())
            self.notes[name].grid_remove()
        self._set_details(())
        self.spectrum.clear()
        self.spectrum_canvas.draw_idle()

    def _make_tabs(self, tables: tuple[TableView, ...]) -> None:
        self.notebook = ttk.Notebook(self.tables_frame)
        self.notebook.pack(fill="both", expand=True)
        for table in tables:
            page = ttk.Frame(self.notebook)
            self.notebook.add(page, text=table.name)
            if table.name == ASSEMBLY:
                self._assembly_page(page)
            else:
                self._table_pane(page, table.name).pack(fill="both", expand=True)

    def _table_pane(self, parent, name: str) -> ttk.Frame:
        """A table with both scroll bars and its note; the tree and note are kept under `name`."""
        pane = ttk.Frame(parent)
        pane.columnconfigure(0, weight=1)
        pane.rowconfigure(0, weight=1)
        tree = ttk.Treeview(pane, show="headings", selectmode="browse", height=5)
        down = ttk.Scrollbar(pane, orient="vertical", command=tree.yview)
        across = ttk.Scrollbar(pane, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=down.set, xscrollcommand=across.set)
        tree.grid(row=0, column=0, sticky="nsew")
        down.grid(row=0, column=1, sticky="ns")
        across.grid(row=1, column=0, sticky="ew")
        note = ttk.Label(pane, style="Muted.TLabel")
        note.grid(row=2, column=0, columnspan=2, sticky="ew")
        note.bind("<Configure>", _wrap_to_width)
        # The tabs exist only once a result is shown, after `bind` has set the presenter.
        tree.bind("<<TreeviewSelect>>", lambda event: self._row_selected(name))
        self.trees[name] = tree
        self.notes[name] = note
        return pane

    def _assembly_page(self, page) -> None:
        """The measure table and the detail lines at the left, the spectrum chart at the right."""
        page.columnconfigure(0, weight=1)
        page.columnconfigure(1, weight=2)
        page.rowconfigure(0, weight=1)
        left = ttk.Frame(page)
        left.grid(row=0, column=0, sticky="nsew")
        left.columnconfigure(0, weight=1)
        left.rowconfigure(0, weight=1)
        self._table_pane(left, ASSEMBLY).grid(row=0, column=0, columnspan=2, sticky="nsew")
        self.details = tkinter.Text(
            left,
            height=4,
            wrap="word",
            state="disabled",
            font="TkDefaultFont",
            background=render.SURFACE,
            foreground=render.INK,
            relief="flat",
            highlightthickness=1,
            highlightbackground=render.RULE,
        )
        scroll = ttk.Scrollbar(left, orient="vertical", command=self.details.yview)
        self.details.configure(yscrollcommand=scroll.set)
        self.details.grid(row=1, column=0, sticky="ew", pady=(4, 0))
        scroll.grid(row=1, column=1, sticky="ns", pady=(4, 0))
        self.spectrum = Figure(facecolor=render.SURFACE)
        self.spectrum_canvas = FigureCanvasTkAgg(self.spectrum, master=page)
        chart = self.spectrum_canvas.get_tk_widget()
        chart.configure(width=240, height=120)  # small, so the chart follows the space it gets
        chart.grid(row=0, column=1, sticky="nsew", padx=(8, 0))

    def _fill_table(self, table: TableView) -> None:
        tree = self.trees[table.name]
        tree.delete(*tree.get_children())
        tree.configure(columns=table.columns)
        font = tkinter.font.nametofont("TkDefaultFont")
        for i, column in enumerate(table.columns):
            cells = (row[i] for row in table.rows[:WIDTH_ROWS])
            tree.heading(column, text=column)
            tree.column(column, width=_column_width(font, column, cells), stretch=False, anchor="w")
        for row in table.rows:
            tree.insert("", "end", values=row)
        note = self.notes[table.name]
        if table.note is None:
            note.grid_remove()
        else:
            note.configure(text=table.note)
            note.grid()

    def _set_details(self, lines) -> None:
        self.details.configure(state="normal")
        self.details.delete("1.0", "end")
        self.details.insert("1.0", "\n".join(lines))
        self.details.configure(state="disabled")

    def _row_selected(self, name: str) -> None:
        """Pass on the row's position: Tk turns numeric-looking cell text into numbers."""
        tree = self.trees[name]
        selection = tree.selection()  # empty while a table is being refilled
        if selection:
            self.presenter.row_selected(name, tree.index(selection[0]))

    def show_busy(self, busy: bool) -> None:
        self.stop_button.state(["!disabled" if busy else "disabled"])
        if busy:
            self.progress.pack(side="right")
            self.progress.start()
        else:
            self.progress.stop()
            self.progress.pack_forget()

    def _graph(self) -> None:
        """The canvas of the plot, in the frame `_explorer` made; empty until a result."""
        self.figure = Figure(facecolor=render.SURFACE)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.canvas_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.configure(width=200, height=150)  # small: it follows the space it gets
        self.canvas_widget.pack(fill="both", expand=True)
        # Not shown: it supplies the zoom and pan modes; exports go through the File menu.
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.canvas_frame, pack_toolbar=False)
        self.left_out_label = ttk.Label(self.strip, style="Muted.TLabel")
        self.left_out_label.pack(side="left", padx=(12, 0))
        self.scene = None
        self.analysis = None
        self._step = 0
        self._home = None  # the whole graph's limits, for Reset view
        self._hovered = None
        self._overview = False
        self._resize_job = None

    def _step_controls(self) -> None:
        self.step_buttons = [ttk.Button(self.steps_frame, text=t) for t in STEP_BUTTONS]
        self.scale = ttk.Scale(self.steps_frame, from_=0, to=0, orient="horizontal")
        self.step_label = ttk.Label(self.steps_frame, width=14)
        first, back, play, forward, last = self.step_buttons
        for widget in (first, back, play, self.scale, forward, last, self.step_label):
            expand = widget is self.scale
            widget.pack(side="left", fill="x" if expand else "none", expand=expand, padx=2)
        self._show_step_label(0, 0)

    def _set_graph_controls(self, enabled: bool) -> None:
        state = "!disabled" if enabled else "disabled"
        for widget in self.view_buttons.values():
            widget.state([state])

    def _show_step_label(self, step: int, total: int) -> None:
        self.step_label.configure(text=f"step {step} / {total}")

    def _measure_aspect(self) -> float:
        """The canvas's width over its height now; sets nothing."""
        self.root.update_idletasks()
        width, height = self.canvas_widget.winfo_width(), self.canvas_widget.winfo_height()
        if width < MIN_CANVAS_PX or height < MIN_CANVAS_PX:
            return render.FIGURE_ASPECT
        return width / height

    def canvas_aspect(self) -> float:
        """Read for a build; the result is laid out again for the canvas as it is then."""
        return self._measure_aspect()

    def show_step(self, step: int, total: int) -> None:
        self._step = step
        self._show_step_label(step, total)
        self.scale.set(step)
        self.scene.show(step)
        self.canvas.draw_idle()

    def set_playing(self, playing: bool) -> None:
        self.step_buttons[2].configure(text="Pause" if playing else "Play")

    def _show_graph(self, result: ResultView) -> None:
        self.analysis = result.analysis
        self.headline.configure(text=result.headline)
        self.scale.configure(to=result.steps)
        if result.layout is None:
            self._build_scene(None)  # past the drawing limit: the canvas holds only the text
        else:
            self._relayout()  # the worker's layout may be for another canvas shape
        for entry in render.legend_entries(result.analysis):
            self._legend_sample(entry)
        self._set_graph_controls(True)

    def _legend_sample(self, entry) -> None:
        sample = ttk.Frame(self.legend_frame)
        sample.pack(side="left", padx=(0, 12))
        line = tkinter.Canvas(
            sample, width=34, height=12, highlightthickness=0, background=render.WINDOW
        )
        dash = DASH_PATTERNS[entry.dash] or None  # None leaves the option out: a solid line
        line.create_line(1, 6, 33, 6, width=2, fill=entry.color, dash=dash)
        line.pack(side="left")
        ttk.Label(sample, text=entry.label).pack(side="left", padx=(4, 0))

    def _build_scene(self, layout) -> None:
        """A new scene at the current step. A new layout starts a new home view."""
        self.scene = render.StepScene(self.figure, self.analysis, layout)
        self.scene.show(self._step)
        axes = self.figure.axes[0]
        self._home = (axes.get_xlim(), axes.get_ylim())
        self.toolbar.update()  # forgets the zoom and pan history of the old axes
        self._overview = len(self.analysis.nodes) > render.OVERVIEW_ABOVE_NODES
        self._hovered = None
        shortened, left_out = self.scene.left_out
        text = f"{shortened} labels shortened, {left_out} left out" if shortened or left_out else ""
        self.left_out_label.configure(text=text)
        self.canvas.draw_idle()

    def _relayout(self) -> None:
        """The one place the window lays out: for the canvas as it is now."""
        self._build_scene(render.compute_layout(self.analysis, self._measure_aspect()))

    def _on_configure(self, event) -> None:
        if self._resize_job is not None:
            self.root.after_cancel(self._resize_job)
        self._resize_job = self.root.after(RESIZE_DELAY_MS, self._resized)

    def _resized(self) -> None:
        self._resize_job = None
        if self.scene is not None:
            self._relayout()

    def _clear_graph(self) -> None:
        self.scene = self.analysis = self._home = None
        self._step = 0
        self.figure.clear()
        self.canvas.draw_idle()
        self.headline.configure(text="")
        self.left_out_label.configure(text="")
        for child in self.legend_frame.winfo_children():
            child.destroy()
        self.scale.configure(to=0)
        self._show_step_label(0, 0)
        self._set_graph_controls(False)

    def _on_motion(self, event) -> None:
        if self.scene is None or not self._overview:
            return
        node = self.scene.node_at(event.x, event.y)
        if node != self._hovered:
            self._hovered = node
            self.scene.show_hover(node)
            self.canvas.draw_idle()

    def _on_scroll(self, event) -> None:
        self._zoom(ZOOM_IN if event.button == "up" else ZOOM_OUT)

    def _zoom(self, factor: float) -> None:
        """Scale the visible range about its centre."""
        if self.scene is None:
            return
        axes = self.figure.axes[0]
        (x0, x1), (y0, y1) = axes.get_xlim(), axes.get_ylim()
        cx, cy, half_w, half_h = (x0 + x1) / 2, (y0 + y1) / 2, (x1 - x0) / 2, (y1 - y0) / 2
        axes.set_xlim(cx - half_w * factor, cx + half_w * factor)
        axes.set_ylim(cy - half_h * factor, cy + half_h * factor)
        self.canvas.draw_idle()

    def _reset_view(self) -> None:
        axes = self.figure.axes[0]
        axes.set_xlim(self._home[0])
        axes.set_ylim(self._home[1])
        self.toolbar.update()
        self.canvas.draw_idle()

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

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
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
ASSEMBLY = "Assembly"
WIDTH_ROWS = 50  # rows of a table read to size its columns
MAX_COLUMN_PX = 320  # a longer cell scrolls instead of widening the table
CELL_PAD_PX = 16
TK_KEY_NAMES = {"Ctrl": "Control", "Esc": "Escape", "Enter": "Return"}


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
        if self.notebook is None:
            self._make_tabs(result.tables)
        for table in result.tables:
            self._fill_table(table)
        self._set_details(result.details)
        render.draw_spectrum(self.spectrum, result.analysis)
        self.spectrum_canvas.draw_idle()

    def clear_result(self) -> None:
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

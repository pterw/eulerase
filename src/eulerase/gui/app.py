"""The window frame: the form, the status bar and the busy state. It decides nothing.

Nothing here runs at import: no widget, no root, no thread. A `Window` is built from a root and
reaches its presenter only through `bind`, so no callback is wired before then.
"""

import tkinter
import tkinter.font
import tkinter.ttk as ttk

from eulerase import render
from eulerase.gui.presenter import MIN_WINDOW, FormValues, fit_window

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

    def show_busy(self, busy: bool) -> None:
        self.stop_button.state(["!disabled" if busy else "disabled"])
        if busy:
            self.progress.pack(side="right")
            self.progress.start()
        else:
            self.progress.stop()
            self.progress.pack_forget()

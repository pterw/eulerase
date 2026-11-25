"""The eulerase command: read the input, analyse it, print the headline and the trails."""

import argparse
import os
import sys
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path

from eulerase import io
from eulerase.analysis import DEFAULT_K, TABLES, Analysis, Settings, analyze, format_cell, rows
from eulerase.errors import EuleraseError
from eulerase.render import FIGURE_FORMATS, save_figure
from eulerase.sequence import Record, parse_kmer_list, parse_records

# Every option that is not an input and acts only on one: (flag, parser destination).
# The outputs are also refused with --gui, whose window writes nothing on the command line.
_OUTPUTS = (
    ("--metrics", "metrics"),
    ("--figure", "figure"),
    ("--csv", "csv"),
    ("--gfa", "gfa"),
    ("--contigs", "contigs"),
    ("--genbank", "genbank"),
)
_NEEDS_INPUT = (
    ("-k", "k"),
    ("--both-strands", "both_strands"),
    ("--unique", "unique"),
    ("--min-count", "min_count"),
    ("--clip-tips", "clip_tips"),
    *_OUTPUTS,
)


class _Parser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        print(f"eulerase: {message}", file=sys.stderr)
        raise SystemExit(2)


def _parser() -> argparse.ArgumentParser:
    p = _Parser(
        prog="eulerase",
        description=(
            "Build the de Bruijn graph of the k-mers of a sequence, a set of reads or a list "
            "of k-mers, and print the Eulerian trails that cover it."
        ),
    )
    source = p.add_mutually_exclusive_group()
    source.add_argument("--seq", nargs="+", metavar="SEQ", help="one or more sequences")
    source.add_argument(
        "--kmers", nargs="+", metavar="KMER", help="a list of k-mers, all one length"
    )
    source.add_argument(
        "--file",
        metavar="PATH",
        help="a FASTA, FASTQ, GenBank or EMBL file, plain or gzip",
    )
    p.add_argument(
        "-k",
        type=int,
        default=None,
        metavar="INT",
        help=f"k-mer length for --seq and --file, at least 2 (default {DEFAULT_K}); "
        "with --kmers, leave -k out: it is taken from the k-mers",
    )
    p.add_argument(
        "--both-strands",
        action="store_true",
        help="also use the reverse complement of every sequence (DNA only; default off)",
    )
    p.add_argument(
        "--unique",
        action="store_true",
        help="treat the input as a set: keep each distinct k-mer once (default off)",
    )
    p.add_argument(
        "--min-count",
        type=int,
        metavar="N",
        help="drop k-mers seen fewer than N times before the contigs are built (default: "
        "keep all). For a set of reads this removes sequencing errors; on one assembled "
        "sequence a low count marks a rarer k-mer, not an error",
    )
    p.add_argument(
        "--clip-tips",
        action="store_true",
        help="before the contigs are built, remove short dead-end branches (under 2k k-mers) "
        "when another branch from the same node is seen more often on average (default off). "
        "For a set of reads this "
        "removes sequencing errors; it can also remove the true end of a sequence. "
        "It removes dead ends, not bubbles: for reads, also give --min-count.",
    )
    p.add_argument("--figure", metavar="PATH", help="write the graph as a png, svg or pdf figure")
    p.add_argument("--csv", metavar="PREFIX", help="write the tables as PREFIX_<table>.csv files")
    p.add_argument("--gfa", metavar="PATH", help="write the contigs as a GFA 1.0 file")
    p.add_argument("--contigs", metavar="PATH", help="write the contigs as a FASTA file")
    p.add_argument("--genbank", metavar="PATH", help="write the contigs as a GenBank file")
    p.add_argument("--metrics", action="store_true", help="print the assembly measures")
    p.add_argument("--gui", action="store_true", help="open the window with the input loaded")
    return p


def _given(args: argparse.Namespace, dest: str) -> bool:
    value = getattr(args, dest)
    return value is not None and value is not False  # 0 is a value the user gave


def _refuse_ignored_options(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """Exit 2 when an option was given that nothing would act on."""
    if args.gui:
        for flag, dest in _OUTPUTS:
            if _given(args, dest):
                parser.error(f"{flag} cannot be used with --gui: the window has its own exports.")
    if args.kmers is not None:
        if args.k is not None:
            parser.error(
                "-k cannot be used with --kmers: leave -k out, the length is taken from the k-mers."
            )
    elif args.seq is None and args.file is None:
        for flag, dest in _NEEDS_INPUT:
            if _given(args, dest):
                parser.error(f"{flag} needs an input: give --seq, --kmers or --file.")


def _records(args: argparse.Namespace) -> list[Record]:
    if args.kmers is not None:
        return parse_kmer_list(" ".join(args.kmers))
    if args.file is not None:
        return io.read_records(args.file)
    return parse_records("\n".join(args.seq))


def _settings(args: argparse.Namespace) -> Settings:
    from_kmers = args.kmers is not None
    return Settings(
        k=args.k if from_kmers or args.k is not None else DEFAULT_K,
        both_strands=args.both_strands,
        unique=args.unique,
        min_multiplicity=args.min_count,
        clip_tips=args.clip_tips,
        input_mode="kmers" if from_kmers else "sequences",
    )


def _form_text(args: argparse.Namespace) -> str:
    """The input as the window's input box holds it."""
    from eulerase.gui.presenter import records_to_text

    if args.kmers is not None:
        return " ".join(args.kmers)
    if args.file is not None:
        return records_to_text(io.read_records(args.file))
    return "\n".join(args.seq)


def _start_window(args: argparse.Namespace, has_input: bool) -> int:
    """Open the window, with the input given loaded into it; the toolkit is imported only here."""
    from eulerase.gui.presenter import FormValues

    form = None
    if has_input:
        form = FormValues(
            text=_form_text(args),
            input_mode="kmers" if args.kmers is not None else "sequences",
            k_text=str(DEFAULT_K if args.k is None else args.k),
            both_strands=args.both_strands,
            unique=args.unique,
            min_multiplicity_text="" if args.min_count is None else str(args.min_count),
            clip_tips=args.clip_tips,
        )
    try:
        from eulerase.gui.app import start
    except ImportError as err:
        if err.name not in ("tkinter", "_tkinter"):
            raise  # some other import failed; the real cause must not be hidden
        raise EuleraseError(
            "The window needs tkinter, which this Python does not have; "
            "the command-line options still work."
        ) from err
    return start(form)


def main(argv: Sequence[str] | None = None) -> int:
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", errors="backslashreplace", newline="\n")
    parser = _parser()
    args = parser.parse_args(argv)
    _refuse_ignored_options(args, parser)
    has_input = any(given is not None for given in (args.seq, args.kmers, args.file))
    try:
        if args.gui or not has_input:
            return _start_window(args, has_input)
        _check_output_paths(args)
        result = analyze(_records(args), _settings(args))
        texts = _export_texts(args, result)
        _print_results(args, result)
        _write_outputs(args, result, texts)
    except EuleraseError as err:
        print(f"eulerase: {err}", file=sys.stderr)
        return 1
    return 0


def _output_paths(args: argparse.Namespace) -> list[tuple[str, str]]:
    """Every file the command would write, as (option, path)."""
    paths = [("--figure", args.figure)] if args.figure is not None else []
    if args.csv is not None:
        paths += [("--csv", path) for path in io.csv_paths(args.csv).values()]
    for option in ("--gfa", "--contigs", "--genbank"):
        if getattr(args, option[2:]) is not None:
            paths.append((option, getattr(args, option[2:])))
    return paths


def _check_output_paths(args: argparse.Namespace) -> None:
    """Refuse the first output that cannot be written, before any work is done."""
    for option, dest in _OUTPUTS:
        if getattr(args, dest) == "":
            raise EuleraseError(f"{option} needs a file name, not an empty path; give one.")
    named: dict[str, str] = {}  # file, as the system sees it -> the option that names it
    for option, path in _output_paths(args):
        file = os.path.normcase(os.path.abspath(path))
        if file in named:
            raise EuleraseError(
                f"{named[file]} and {option} name the same file {path}; "
                "give each output its own path."
            )
        named[file] = option
        target = Path(path)
        folder = target.parent
        if option == "--figure" and target.suffix.lower().lstrip(".") not in FIGURE_FORMATS:
            formats = ", ".join(f".{f}" for f in FIGURE_FORMATS)
            raise EuleraseError(
                f"--figure '{path}' is not a figure format; end the name in {formats}."
            )
        if target.is_dir():
            raise EuleraseError(f"{option} '{path}' is a folder; give a file name instead.")
        if not folder.is_dir():
            raise EuleraseError(
                f"{option} '{path}' is in a folder that does not exist; create it first."
            )
        overwrite = target.exists()
        if not os.access(target if overwrite else folder, os.W_OK):
            what = (
                "is a file that cannot be written"
                if overwrite
                else "is in a folder that cannot be written"
            )
            raise EuleraseError(f"{option} '{path}' {what}; choose another.")


@contextmanager
def _writing(path: str) -> Iterator[None]:
    """Report an OSError raised while writing `path` as one sentence."""
    try:
        yield
    except OSError as err:
        raise EuleraseError(f"cannot write {path}: {err.strerror or err}") from err


def _print_results(args: argparse.Namespace, result: Analysis) -> None:
    print(result.headline)
    for note in result.notes:
        print(f"Note: {note}")
    for i, (trail, spelled) in enumerate(zip(result.trails, result.spelled, strict=True), 1):
        print(f"trail {i} ({len(trail)} edges): {spelled}")
    if args.metrics:
        _print_metrics(result)


def _export_texts(args: argparse.Namespace, result: Analysis) -> dict[str, str]:
    """Build the GFA, FASTA and GenBank texts asked for; a refusal comes before any output."""
    writers = {"gfa": io.gfa_text, "contigs": io.fasta_text, "genbank": io.genbank_text}
    return {
        dest: write(result) for dest, write in writers.items() if getattr(args, dest) is not None
    }


def _write_outputs(args: argparse.Namespace, result: Analysis, texts: dict[str, str]) -> None:
    """Write each requested file in the order figure, csv, gfa, contigs, genbank."""
    if args.figure is not None:
        with _writing(args.figure):
            save_figure(result, args.figure)
        print(f"Wrote {args.figure}")
    if args.csv is not None:
        paths = io.csv_paths(args.csv)
        for table in TABLES:
            path = paths[table.lower()]
            with _writing(path):
                io.write_text_file(path, io.csv_text(result, table))
            print(f"Wrote {path}")
    for dest, text in texts.items():
        path = getattr(args, dest)
        with _writing(path):
            io.write_text_file(path, text)
        print(f"Wrote {path}")


def _print_metrics(result: Analysis) -> None:
    for measure, value in rows(result, "Assembly"):
        print(f"{measure}: {value}")
    for contig in result.contigs:
        print(
            f"contig {contig.id} length {contig.length} "
            f"multiplicity {format_cell(contig.mean_multiplicity)} "
            f"cycle {format_cell(contig.circular)}"
        )

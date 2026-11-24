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

_NOT_BUILT = "The window is not built yet; give --seq, --kmers or --file."

# Every option that is not an input and acts only on one: (flag, parser destination).
_NEEDS_INPUT = (
    ("-k", "k"),
    ("--both-strands", "both_strands"),
    ("--unique", "unique"),
    ("--min-count", "min_count"),
    ("--clip-tips", "clip_tips"),
    ("--metrics", "metrics"),
    ("--figure", "figure"),
    ("--csv", "csv"),
    ("--gfa", "gfa"),
    ("--contigs", "contigs"),
    ("--genbank", "genbank"),
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
    p.add_argument("--gui", action="store_true", help="open the window (not built yet)")
    return p


def _refuse_ignored_options(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """Exit 2 when an option was given that nothing would act on."""
    if args.kmers is not None:
        if args.k is not None:
            parser.error(
                "-k cannot be used with --kmers: leave -k out, the length is taken from the k-mers."
            )
    elif args.seq is None and args.file is None:
        for flag, dest in _NEEDS_INPUT:
            value = getattr(args, dest)
            if value is not None and value is not False:  # 0 is a value the user gave
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


def main(argv: Sequence[str] | None = None) -> int:
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    parser = _parser()
    args = parser.parse_args(argv)
    _refuse_ignored_options(args, parser)
    if args.gui or all(given is None for given in (args.seq, args.kmers, args.file)):
        print(_NOT_BUILT, file=sys.stderr)
        return 1
    try:
        _check_output_paths(args)
        result = analyze(_records(args), _settings(args))
        _print_results(args, result)
        _write_outputs(args, result)
    except EuleraseError as err:
        print(f"eulerase: {err}", file=sys.stderr)
        return 1
    return 0


def _output_paths(args: argparse.Namespace) -> list[tuple[str, str]]:
    """Every file the command would write, as (option, path)."""
    paths = [("--figure", args.figure)] if args.figure else []
    if args.csv:
        paths += [("--csv", path) for path in io.csv_paths(args.csv).values()]
    for option in ("--gfa", "--contigs", "--genbank"):
        if getattr(args, option[2:]):
            paths.append((option, getattr(args, option[2:])))
    return paths


def _check_output_paths(args: argparse.Namespace) -> None:
    """Refuse the first output that cannot be written, before any work is done."""
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
        if not os.access(folder, os.W_OK):
            raise EuleraseError(
                f"{option} '{path}' is in a folder that cannot be written; choose another."
            )


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


def _write_outputs(args: argparse.Namespace, result: Analysis) -> None:
    """Write each requested file in the order figure, csv, gfa, contigs, genbank."""
    if args.figure:
        with _writing(args.figure):
            save_figure(result, args.figure)
        print(f"Wrote {args.figure}")
    if args.csv:
        paths = io.csv_paths(args.csv)
        for table in TABLES:
            path = paths[table.lower()]
            with _writing(path):
                io.write_text_file(path, io.csv_text(result, table))
            print(f"Wrote {path}")
    for path, writer in (
        (args.gfa, io.gfa_text),
        (args.contigs, io.fasta_text),
        (args.genbank, io.genbank_text),
    ):
        if path:
            with _writing(path):
                io.write_text_file(path, writer(result))
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

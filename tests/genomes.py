"""Test support: random genomes, an independent reverse complement, brute-force k-mers."""

import os
import random
from collections.abc import Sequence
from pathlib import Path

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def random_dna(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s: str) -> str:
    return s.translate(_COMPLEMENT)[::-1]


def least_rotation(s: str) -> str:
    """The lexicographically least rotation of a circle, by trying every rotation."""
    return min(s[i:] + s[:i] for i in range(len(s)))


def is_the_genome(sequences: Sequence[str], genome: str, both: bool) -> bool:
    """True when the sequences are the genome alone (or, both strands, its reverse complement)."""
    allowed = {genome, revcomp(genome)} if both else {genome}
    return len(sequences) == 1 and sequences[0] in allowed


def brute_kmers(sequences: Sequence[str], k: int, both: bool) -> list[str]:
    strands = [s for seq in sequences for s in ((seq, revcomp(seq)) if both else (seq,))]
    return [s[i : i + k] for s in strands for i in range(len(s) - k + 1)]


def child_env() -> dict[str, str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1] / "src")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return env

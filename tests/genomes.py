"""Test support: random genomes, an independent reverse complement, brute-force k-mers."""

import os
import random
from collections.abc import Sequence
from pathlib import Path

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def random_dna(n: int, seed: int) -> str:
    return "".join(random.Random(seed).choice("ACGT") for _ in range(n))


def revcomp(s: str) -> str:
    return s.translate(_COMPLEMENT)[::-1]


def brute_kmers(sequences: Sequence[str], k: int, both: bool) -> list[str]:
    strands = [s for seq in sequences for s in ((seq, revcomp(seq)) if both else (seq,))]
    return [s[i : i + k] for s in strands for i in range(len(s) - k + 1)]


def child_env() -> dict[str, str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1] / "src")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return env

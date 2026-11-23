"""Third gate: run the tests on a copy of the tracked working-tree files, and nothing else.

A test that passes only because of an untracked file, an installed copy of the package or
the original source tree fails here. Usage: python tools/check_fresh_export.py [REPO]
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def tracked_files(repo: Path) -> list[Path]:
    listing = subprocess.run(
        ["git", "ls-files", "-z"], cwd=repo, check=True, capture_output=True, timeout=60
    ).stdout
    return [Path(os.fsdecode(name)) for name in listing.split(b"\0") if name]


def copy_tracked(repo: Path, target: Path) -> None:
    for relative in tracked_files(repo):
        source = repo / relative
        if source.is_file():
            (target / relative).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target / relative)


def fresh_environment(copy: Path) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "EULERASE_SOURCE")}
    env["PYTHONPATH"] = str(copy / "src")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return env


def package_is_in_copy(copy: Path, env: dict[str, str]) -> bool:
    probe = subprocess.run(
        [sys.executable, "-c", "import eulerase; print(eulerase.__file__)"],
        cwd=copy,
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
    )
    return probe.returncode == 0 and Path(probe.stdout.strip()).resolve().is_relative_to(
        copy.resolve()
    )


def main(repo: str = ".") -> int:
    with tempfile.TemporaryDirectory() as scratch:
        copy = Path(scratch)
        copy_tracked(Path(repo), copy)
        env = fresh_environment(copy)
        if not package_is_in_copy(copy, env):
            print("eulerase did not import from the export: is src/eulerase tracked?")
            return 1
        tests = subprocess.run(
            [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"],
            cwd=copy,
            env=env,
            timeout=1800,
        )
        return tests.returncode


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:2]))

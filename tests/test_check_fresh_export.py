"""The fresh-export gate runs the tests on a copy of the tracked working-tree files."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

TOOL = Path(__file__).resolve().parents[1] / "tools" / "check_fresh_export.py"

pytestmark = pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")

PASSING = "def test_ok():\n    assert True\n"
FAILING = "def test_ok():\n    assert False\n"
NEEDS_NEW_FILE = (
    "from pathlib import Path\n\n\n"
    "def test_ok():\n"
    "    assert (Path(__file__).parent / 'new.txt').exists()\n"
)


def git_env() -> dict[str, str]:
    env = dict(os.environ)
    for name in ("GIT_AUTHOR", "GIT_COMMITTER"):
        env[f"{name}_NAME"] = "Test"
        env[f"{name}_EMAIL"] = "nobody"
    env["GIT_CONFIG_GLOBAL"] = os.devnull
    env["GIT_CONFIG_NOSYSTEM"] = "1"
    return env


def git(repo: Path, *args: str) -> None:
    subprocess.run(["git", *args], cwd=repo, env=git_env(), check=True, capture_output=True)


def make_repo(root: Path, test_source: str, *, with_package: bool = True) -> Path:
    repo = root / "repo"
    (repo / "tests").mkdir(parents=True)
    (repo / "tests" / "test_x.py").write_text(test_source, encoding="utf-8")
    if with_package:
        (repo / "src" / "eulerase").mkdir(parents=True)
        (repo / "src" / "eulerase" / "__init__.py").write_text("", encoding="utf-8")
    git(repo, "init", "-q")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "initial")
    return repo


def run_gate(repo: Path, **extra_env: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(TOOL), str(repo)],
        capture_output=True,
        text=True,
        timeout=120,
        env={**git_env(), **extra_env},
    )


def verdict(result: subprocess.CompletedProcess) -> str:
    """pytest's own summary line, or the exit status when pytest never ran."""
    lines = result.stdout.strip().splitlines()
    summary = lines[-1] if lines else ""
    return (
        "passed"
        if result.returncode == 0 and "1 passed" in summary
        else ("failed" if result.returncode != 0 and "1 failed" in summary else f"other: {result}")
    )


def test_passing_repository_exits_zero(tmp_path):
    assert verdict(run_gate(make_repo(tmp_path, PASSING))) == "passed"


@pytest.mark.parametrize("outside", ["untracked-file", "module-on-inherited-pythonpath"])
def test_a_test_that_needs_something_outside_the_export_fails(tmp_path, outside):
    if outside == "untracked-file":
        repo = make_repo(tmp_path, NEEDS_NEW_FILE.replace("new.txt", "untracked.txt"))
        (repo / "tests" / "untracked.txt").write_text("x", encoding="utf-8")
        assert verdict(run_gate(repo)) == "failed"
    else:
        repo = make_repo(tmp_path, "def test_ok():\n    import outside_helper\n")
        elsewhere = tmp_path / "elsewhere"
        elsewhere.mkdir()
        (elsewhere / "outside_helper.py").write_text("", encoding="utf-8")
        assert verdict(run_gate(repo, PYTHONPATH=str(elsewhere))) == "failed"


def test_repository_without_the_package_fails_before_the_tests_run(tmp_path):
    result = run_gate(make_repo(tmp_path, PASSING, with_package=False))
    assert result.returncode != 0
    assert "passed" not in result.stdout


@pytest.mark.parametrize(
    ("committed", "working_tree", "expected"),
    [(PASSING, FAILING, "failed"), (FAILING, PASSING, "passed")],
    ids=["uncommitted-edit-breaks", "uncommitted-edit-fixes"],
)
def test_gate_reads_the_working_tree_not_the_commit(tmp_path, committed, working_tree, expected):
    repo = make_repo(tmp_path, committed)
    (repo / "tests" / "test_x.py").write_text(working_tree, encoding="utf-8")
    assert verdict(run_gate(repo)) == expected


@pytest.mark.parametrize(("staged", "expected"), [(True, "passed"), (False, "failed")])
def test_only_a_staged_new_file_is_in_the_export(tmp_path, staged, expected):
    repo = make_repo(tmp_path, NEEDS_NEW_FILE)
    (repo / "tests" / "new.txt").write_text("x", encoding="utf-8")
    if staged:
        git(repo, "add", "tests/new.txt")
    assert verdict(run_gate(repo)) == expected

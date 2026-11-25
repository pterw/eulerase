"""The gates workflow runs the gates the repository documents.

These tests compare the workflow with the gate commands in AGENTS.md and with
the shape the workflow is held to. They do not check what the gates do.
"""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
INSTALL = 'python -m pip install -e ".[dev]"'

FROZEN_TRIGGER = [
    "on:",
    "  pull_request:",
    "    paths-ignore:",
    "      - docs/**/*.md",
    "      - CONTRIBUTING.md",
    "      - .github/ISSUE_TEMPLATE/*.md",
    "      - .github/PULL_REQUEST_TEMPLATE.md",
    "  push:",
    "    branches: [main]",
]
FROZEN_PERMISSIONS = ["permissions:", "  contents: read"]
TOP_KEYS = {"name", "on", "permissions", "concurrency", "jobs"}
JOB_KEYS = {"name", "runs-on", "strategy", "timeout-minutes", "steps"}
STEP_KEYS = [
    ["uses", "with"],
    ["uses", "with"],
    ["name", "run"],
    ["name", "run"],
    ["name", "run"],
]
JOB_NAME = "    name: gates (${{ matrix.os }}, ${{ matrix.python-version }})"
PYTHON_FROM_MATRIX = "          python-version: ${{ matrix.python-version }}"
FROZEN_MATRIX = [
    "        os: [windows-latest, ubuntu-latest, macos-latest]",
    '        python-version: ["3.13"]',
    "        include:",
    "          - os: ubuntu-latest",
    '            python-version: "3.15"',
]
WITH_KEYS = [
    ["persist-credentials"],
    ["cache", "cache-dependency-path", "python-version"],
]


def documented_gates(agents_text: str) -> list[str]:
    """The `run:` values of the fenced `gates` block's items not marked `local: yes`."""
    lines = [line.strip() for line in agents_text.splitlines()]
    assert "```gates" in lines, "AGENTS.md has no gates block"
    block = lines[lines.index("```gates") + 1 :]
    block = block[: block.index("```")]
    items: list[list[str]] = []
    for line in block:
        if line.startswith("- name:"):
            items.append([])
        items[-1].append(line)
    return [
        next(x[4:].strip() for x in item if x.startswith("run:"))
        for item in items
        if "local: yes" not in item
    ]


def workflow_runs(workflow_text: str) -> list[str]:
    """The value of every `run:` line; a block or empty value is refused."""
    runs = []
    for line in workflow_text.splitlines():
        line = line.strip().removeprefix("- ")
        if line.startswith("run:"):
            value = line[4:].strip()
            assert value not in ("", "|", ">"), f"run is not one line: {line}"
            runs.append(value)
    return runs


def _clean(text: str) -> list[str]:
    """Lines without comments and blanks, keeping indentation."""
    lines = [line.split(" #")[0].rstrip() for line in text.splitlines()]
    return [x for x in lines if x.strip() and not x.lstrip().startswith("#")]


def _indent(line: str) -> int:
    return len(line) - len(line.lstrip())


def _block(lines: list[str], header: str) -> list[str]:
    """The lines under `header`, which are indented deeper than it."""
    if header not in lines:
        return []
    rest = lines[lines.index(header) + 1 :]
    deeper = []
    for line in rest:
        if _indent(line) <= _indent(header):
            break
        deeper.append(line)
    return deeper


def _keys(lines: list[str], indent: int) -> list[str]:
    return [x.strip().split(":")[0] for x in lines if _indent(x) == indent]


def _steps(job: list[str]) -> list[list[str]]:
    steps: list[list[str]] = []
    for line in _block(job, "    steps:"):
        if line.startswith("      - "):
            steps.append([])
            line = "        " + line[8:]
        steps[-1].append(line)
    return steps


def workflow_flaws(workflow_text: str) -> list[str]:
    """Every way the workflow leaves the shape it is held to."""
    flaws = []
    raw = [line.rstrip() for line in workflow_text.splitlines()]
    start = raw.index("on:") if "on:" in raw else len(raw)
    end = next((i for i in range(start, len(raw)) if raw[i] in ("permissions:", "jobs:")), len(raw))
    if raw[start:end] != FROZEN_TRIGGER:
        flaws.append("trigger differs from the frozen block")

    lines = _clean(workflow_text)
    if set(_keys(lines, 0)) - TOP_KEYS:
        flaws.append("unexpected top-level key")
    if ["permissions:", *_block(lines, "permissions:")] != FROZEN_PERMISSIONS:
        flaws.append("permissions are not contents: read")
    if any("{" in re.sub(r"\$\{\{.*?\}\}", "", x) for x in lines):
        flaws.append("flow-style braces")

    jobs = _block(lines, "jobs:")
    if _keys(jobs, 2) != ["gates"]:
        flaws.append("jobs are not exactly gates")
    job = _block(jobs, "  gates:")
    if not set(_keys(job, 4)) <= JOB_KEYS:
        flaws.append("unexpected job key")
    if "    runs-on: ${{ matrix.os }}" not in job:
        flaws.append("runs-on is not the matrix system")

    strategy = _block(job, "    strategy:")
    if not set(_keys(strategy, 6)) <= {"fail-fast", "matrix"}:
        flaws.append("unexpected strategy key")
    if "      fail-fast: false" not in strategy:
        flaws.append("fail-fast is not false")
    matrix = _block(strategy, "      matrix:")
    if matrix != FROZEN_MATRIX:
        flaws.append("matrix is not the four frozen jobs")
    if JOB_NAME not in job:
        flaws.append("job name does not tell the four jobs apart")

    steps = _steps(job)
    if [sorted(_keys(s, 8)) for s in steps] != STEP_KEYS:
        flaws.append("steps are not the five frozen steps")
    elif [sorted(_keys(_block(s, "        with:"), 10)) for s in steps[:2]] != WITH_KEYS:
        flaws.append("a uses step has other with keys")
    elif PYTHON_FROM_MATRIX not in _block(steps[1], "        with:"):
        flaws.append("setup-python does not take the matrix's Python")
    if any(len(s) != 2 for s in steps if "run" in _keys(s, 8)):
        flaws.append("a gate step has a continued or extra line")
    return flaws


def problems(workflow_text: str, documented: list[str]) -> list[str]:
    """Every difference between the workflow and the documented gates."""
    try:
        runs = workflow_runs(workflow_text)
    except AssertionError as error:
        return [str(error)]
    found = []
    if runs[:1] != [INSTALL]:
        found.append("install line differs")
    if runs[1:] != documented:
        found.append("gate commands differ from the documented ones")
    return found + workflow_flaws(workflow_text)


def _read(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def test_workflow_runs_the_documented_gates_in_order():
    documented = documented_gates(_read("AGENTS.md"))
    assert len(documented) >= 2
    assert problems(_read(".github/workflows/gates.yml"), documented) == []


def _edit(old: str, new: str):
    def apply(text: str) -> str:
        assert old in text, f"the frozen workflow has no {old!r}"
        return text.replace(old, new, 1)

    return apply


def _swap_first_two_gates(text: str) -> str:
    first, second = "python -m pytest -q -p no:cacheprovider", "python -m ruff check ."
    assert first != second
    return text.replace(first, "@@").replace(second, first).replace("@@", second)


RUFF_STEP = "      - name: ruff\n        run: python -m ruff check .\n"
EXPORT_STEP = "      - name: fresh-export\n        run: python tools/check_fresh_export.py\n"
MACOS_LIST = "[windows-latest, ubuntu-latest, macos-latest]"
INCLUDE = '        include:\n          - os: ubuntu-latest\n            python-version: "3.15"\n'
DRIFTS = {
    "gate command altered": _edit("python -m ruff check .", "python -m ruff check src"),
    "gate not marked local missing": _edit(RUFF_STEP, ""),
    "fresh-export step present": _edit(RUFF_STEP, RUFF_STEP + EXPORT_STEP),
    "macOS entry removed": _edit(MACOS_LIST, "[windows-latest, ubuntu-latest]"),
    "macOS left only in a comment": _edit(
        MACOS_LIST, "[windows-latest, ubuntu-latest]  # macos-latest"
    ),
    "3.15 include removed": _edit(INCLUDE, ""),
    "3.15 include moved to windows": _edit(
        "- os: ubuntu-latest\n            python-version",
        "- os: windows-latest\n            python-version",
    ),
    "python-version hard-coded in setup-python": _edit(
        PYTHON_FROM_MATRIX, '          python-version: "3.13"'
    ),
    "gates swapped": _swap_first_two_gates,
    "extra run between gates": _edit(
        RUFF_STEP, "      - name: x\n        run: echo x\n" + RUFF_STEP
    ),
    "run block": _edit(
        "        run: python -m ruff check .\n",
        "        run: |\n          python -m ruff check .\n",
    ),
    "seventh step before the gates": _edit(
        "      - name: pytest\n",
        "      - name: x\n        run: python -m pip list\n      - name: pytest\n",
    ),
    "second command on the install line": _edit(INSTALL, INSTALL + "; echo hidden"),
    "src added to paths-ignore": _edit(
        "      - docs/**/*.md\n", "      - docs/**/*.md\n      - src/**\n"
    ),
    "docs directory ignored": _edit("      - docs/**/*.md\n", "      - docs/**\n"),
    "branches-ignore on push": _edit("    branches: [main]", "    branches-ignore: [main]"),
    "push trigger removed": _edit("  push:\n    branches: [main]\n", ""),
    "pull_request types": _edit("  pull_request:\n", "  pull_request:\n    types: [closed]\n"),
    "permissions widened": _edit("contents: read", "contents: write"),
    "exclude entry": _edit(INCLUDE, INCLUDE + "        exclude:\n          - os: macos-latest\n"),
    "hard-coded runs-on": _edit("runs-on: ${{ matrix.os }}", "runs-on: ubuntu-latest"),
    "fail-fast true": _edit("fail-fast: false", "fail-fast: true"),
    "continue-on-error on a gate": _edit(
        RUFF_STEP, RUFF_STEP + "        continue-on-error: true\n"
    ),
    "working-directory on a gate": _edit(
        RUFF_STEP, RUFF_STEP + "        working-directory: docs\n"
    ),
    "if on the job": _edit(
        "    runs-on: ${{ matrix.os }}\n", "    runs-on: ${{ matrix.os }}\n    if: false\n"
    ),
    "continued run line": _edit(
        "        run: python -m pytest -q -p no:cacheprovider\n",
        "        run: python -m pytest -q -p no:cacheprovider\n          || true\n",
    ),
    "flow-style with": _edit(
        "        with:\n          persist-credentials: false\n",
        "        with: {persist-credentials: false}\n",
    ),
}


@pytest.mark.parametrize("drift", DRIFTS.values(), ids=DRIFTS.keys())
def test_oracle_rejects_a_drifted_workflow(drift):
    frozen = _read(".github/workflows/gates.yml")
    documented = documented_gates(_read("AGENTS.md"))
    assert problems(frozen, documented) == []
    broken = drift(frozen)
    assert broken != frozen
    assert problems(broken, documented) != []


def test_the_matrix_drift_cases_differ_from_each_other_and_from_the_real_text():
    frozen = _read(".github/workflows/gates.yml")
    macos = DRIFTS["macOS entry removed"](frozen)
    python = DRIFTS["3.15 include removed"](frozen)
    assert len({frozen, macos, python}) == 3

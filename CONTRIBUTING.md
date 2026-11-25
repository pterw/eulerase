# Contributing to Eulerase

Eulerase is a desktop and command-line tool for learning and exploring de Bruijn graphs: Eulerian
trails and contigs from short sequences. The design, and what is in and out of scope, is in
[the design spec](docs/superpowers/specs/2026-11-09-eulerase-v0.1-design.md).

## Setup

Python 3.13. From a clone, in a virtual environment:

```
python -m venv .venv
```

Activate it (`.venv\Scripts\activate` on Windows, `source .venv/bin/activate` elsewhere), then
install:

```
python -m pip install -e ".[dev]"
```

The quotes around the extras matter: they keep the shell from treating `[dev]` as a pattern.

## The gates

Three gates run in this order from the repository root, before every commit:

```
python -m pytest -q -p no:cacheprovider
python -m ruff check .
python tools/check_fresh_export.py
```

The third runs the tests on a copy of the tracked working-tree files only, so a new file must be
staged (`git add <file>`) before it runs, or the copy lacks it.

CI runs the same three gates on Windows and Linux for every push to `main` and for every pull
request that changes anything but documentation. The documentation paths are the four in the
workflow's `paths-ignore` (Markdown under `docs/`, `CONTRIBUTING.md`, the issue templates and the
pull request template); they hold no code, test, tool, `pyproject.toml`, `README.md`, `LICENSE`,
`AGENTS.md` or workflow. A documentation-only pull request therefore shows no CI run.

## Releases

The version lives in `pyproject.toml` and in `src/eulerase/__init__.py`, and a test holds the two
equal. A release is a pull request that changes the version and adds its entry to
[the changelog](CHANGELOG.md).

After the pull request is merged and the gates are green on `main`, the release workflow tags the
commit those gates tested `v<version>` and attaches the source distribution and the wheel to a
GitHub release of that name. A merge that does not change the version releases nothing. After a
release pull request is merged, wait for its run to finish before merging anything else.

Nothing is uploaded to a package index, and no one tags by hand.

## Commits

Use `type(scope): summary`, for example `fix(orfs): ...`. One concern per commit. Stage files by
name. The body is for a developer reading `git log`: the problem with a concrete case, the change,
why this design, and which tests pin it. Put no task, finding or review identifiers in a message,
and no attribution trailer. Never amend, rebase, force-push or skip the hooks.

## Tests

The rules, from section 8.1 of the spec:

- A test names the behaviour it guards and fails if that behaviour breaks. It checks an intended
  result against an independent answer (brute force, a direct library call, a mathematical
  identity), never the code against itself.
- One test per behaviour. Near-duplicates become one parametrized test. No test pins wording or a
  constant.
- No test needs a display, the network or the home directory. No test opens a window.
- A wrong green is worse than a red.

## Changing behaviour

The spec is the authority. A change to behaviour starts there: where the code and the spec
disagree, the spec wins until it is changed. A departure is made by changing the spec in a commit of
its own, with the reason. A parked feature (section 3 of the spec) returns only by passing its truth
tests (section 8.3), never because it is easy to build.

## Pull requests

Fill in the [pull request template](.github/PULL_REQUEST_TEMPLATE.md): what changed and why, the
head SHA, each gate's last output line quoted (not paraphrased), and under "Not run" anything you did
not run or check. To report a problem or ask for a feature, use the
[bug report](.github/ISSUE_TEMPLATE/bug_report.md) or
[feature request](.github/ISSUE_TEMPLATE/feature_request.md) template.

## Agents and models

The rules an agent or model follows when it works on the repository are in [AGENTS.md](AGENTS.md).

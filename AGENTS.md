# AGENTS.md

Engineering rules for anyone, human or model, who works on Eulerase. The design it follows is
[the design spec](docs/superpowers/specs/2026-11-09-eulerase-v0.1-design.md). A section number
below, and the words "this document" where a rule is copied from the spec, refer to that spec.

## Document roles

| Document | What it is for |
|---|---|
| [the design spec](docs/superpowers/specs/2026-11-09-eulerase-v0.1-design.md) | What the tool is, what it does and does not do, the defects to fix, the tests it needs. It is the authority: where another document disagrees with it, the spec wins until it is changed. |
| [the plans](docs/superpowers/plans) | How the work is carried out, step by step. A plan serves the spec and never overrides it. |
| [README.md](README.md) | What a user needs: what the tool is, how to install and run it, the limits, what it does not do. It describes only what is shipped. |
| [CONTRIBUTING.md](CONTRIBUTING.md) | How a person makes a change: setup, the gates, commits, pull requests. |
| [CHANGELOG.md](CHANGELOG.md) | The changes of each release. |
| [SECURITY.md](SECURITY.md) | How to report a vulnerability. |
| AGENTS.md | The rules an agent or model follows when it works on the repository. |

Where a document describes what the code does and the code does something else, the code is
checked and the document is corrected, unless the spec says the code is wrong.

## Setup

Python 3.13. From a clone, in a virtual environment:

```
python -m pip install -e ".[dev]"
```

The quotes around the extras keep the shell from treating `[dev]` as a pattern.

## Gates

Three gates run in this order from the repository root, before every commit. A gate passes when
its command exits 0, and a report quotes the last line of its output. This block is the
repository's statement of the gate commands. The workflow runs every gate that is not marked
`local: yes`, and its commands must equal the block's.

```gates
- name: pytest
  run: python -m pytest -q -p no:cacheprovider
  pass: exit 0
  quote: last line
- name: ruff
  run: python -m ruff check .
  pass: exit 0
  quote: last line
- name: fresh-export
  run: python tools/check_fresh_export.py
  pass: exit 0
  quote: last line
  local: yes
```

The third gate runs the tests on a copy of the tracked working-tree files, so a new file must be
staged by name before it runs, or the copy lacks it. The core imports neither the window toolkit
nor the plotting library, and nothing imports pandas, scikit-learn or scipy; a test and the linter
enforce both.

## Commits

`type(scope): summary`, one concern per commit. The body is for a developer reading `git log`: the
problem with a concrete case, the change, why this design, and which tests pin it. Put no task,
finding or review identifier in a message, and no attribution trailer. Stage files by name. Never
amend, rebase, force-push or skip the hooks.

## Releases

The version is in `pyproject.toml` and in `src/eulerase/__init__.py`; a test holds the two equal.
A release is a pull request that changes the version and adds its entry to `CHANGELOG.md`. After
it is merged and the gates are green on `main`, the release workflow tags the commit its gates
tested `v<version>` and attaches a source distribution and a wheel to a GitHub release. A merge
that does not change the version releases nothing. No one tags by hand. The release workflow also
publishes each release, after its GitHub release, to TestPyPI and then PyPI.

## Markdown

- Write plain English in short sentences. Headings are in sentence case.
- Write dates as ISO 8601 (`2026-11-24`).
- A document holds no task, finding or review identifier. Write an anti-pattern in words.
- Links between files in the repository are relative. The only absolute addresses are the
  repository's own address (the package metadata, the security policy's link for a private report
  and the status badge at the top of the README), the three static badges of that row, and the SVG
  namespace of the two logo files.

## When you change, also update

| When you change | Also update |
|---|---|
| The gate commands | The `## Gates` block above, `CONTRIBUTING.md` and `.github/workflows/gates.yml` (a test compares the workflow with the block and with the shape the workflow must keep). |
| A limit in the code | The README's table of limits. |
| A command-line option or its help text | The README's command-line section. |
| The palette in `eulerase.render` | The two logo files in `docs/assets`, which use it. |
| The version | `pyproject.toml` and `src/eulerase/__init__.py` (a test holds them equal), and the `CHANGELOG.md` entry. |

## Work model

- The owner merges and pushes.
- A change starts from a written task or issue, and the spec governs it.
- An agent never pushes, amends, rebases or force-pushes.

## Tests

The rules of section 8.1 of the spec:

- A test names the behaviour it guards and fails if that behaviour breaks. It checks an intended
  result against an independent answer, never the code against itself.
- One test per behaviour. Near-duplicates become one parametrized test. No test pins wording or a
  constant.
- No test needs a display, the network, or the user's home directory.
- A wrong green is worse than a red.

No test opens a window.

## The window

Each statement below is a rule of the spec or was read from the source.

- The presenter, `eulerase.gui.presenter`, holds every decision and imports neither tkinter nor
  matplotlib. `eulerase.gui.app` is the view and nothing else.
- No test opens a window or constructs `tk.Tk()`: `tests/conftest.py` replaces `tkinter.Tk` and
  `tkinter.Toplevel` with classes that raise before the toolkit runs. A test that needs tkinter
  skips with a stated reason where it is missing.
- The hands-on pass is a person's, not a test.
- The palette is defined once, in `eulerase.render`. No colour literal appears elsewhere in the
  window code.
- The current position, the starts and the ends never rely on colour alone: an outline, a ring or
  a word marks them. Trails 1 to 3 differ by colours that pass the colour-vision check, with the
  legend, the trail numbers and the Edges table beside them. From trail 4 on a dash pattern is
  added.
- The layout and the placement of labels use no random number and no font metrics, so the same
  input and canvas shape give the same coordinates.
- Only the window's own thread touches a widget. A build runs in a worker that returns its outcome
  through a queue the window polls; a running build can be abandoned, and its result is then
  discarded.

## Anti-patterns

Do none of these.

- A default that removes data (an automatic cutoff).
- A filter that compares nothing and depends on input order.
- Ordering taken from hashing or from a set.
- A label or option name that says "coverage" or "depth" for a count of occurrences.
- A check that reads the code against itself.
- A test that pins wording or a constant.
- Output taken from the clock, the host or a path.
- Text files written with the platform's line endings.
- A format guessed from a file name.
- A path that does not exist read as a sequence.
- An option accepted and ignored.
- A failure that ends in a traceback after some output was written.
- An importing side effect.
- A table whose header and value come from two definitions.
- A failed step that leaves an old result visible.
- A count reported as depth, copy number or evidence of anything in an organism.
- A parked feature built because it was easy.
- A decision made in the window instead of the presenter.
- A test that opens a window.
- A colour literal outside the palette constants.
- A widget touched from a worker thread.
- A position on a figure that depends on a random number or on font metrics.

## Picking up the work

The rules of section 12.1 of the spec, for any agent or model that picks the work up:

- **One task, one commit**, or the commits its steps name. Never gather several tasks into one
  commit: a change that cannot be rolled back or reviewed alone is a defect of process.
- **Lean tests** (section 8.1). Write the tests a task names. Add one only where an intended
  behaviour would otherwise have no test, and say so. More tests are not more safety: a test that
  pins wording, repeats another, or checks the code against itself is removed.
- **An unclear instruction is not a conflict.** Resolve it from the section of this document the
  task cites, say in the task's report which reading you took and why, and carry on. Stop and ask
  only for a real contradiction between two rules, or where the choice would change a number, a
  label or a file format a user sees.
- **Carry, then repair.** The verified code is moved as it is; each defect is repaired in its own
  step with its failing test first. Do not rewrite what works.
- **Say what was run.** A report lists the commands that were run with their output, and
  separately what was not run.

## Invariants

1. Every number shown has a test against an independent oracle.
2. Nothing is removed from the user's data unless the user asked, and the result says what was
   removed.
3. A label says what was computed, in plain words.
4. The same input gives the same output.
5. A failure is a message. A failed step leaves no stale result behind.
6. The core imports neither the GUI nor the plotting library.
7. k-mers from one sequence are not independent observations. No statistic treats them as
   replicates, and no result is shown for a design that cannot support it.
8. A multiplicity is a count of occurrences in the input. It is never shown as depth, copy number or
   evidence of anything in an organism (section 5.1).
9. A parked feature returns only by passing its truth tests (section 8.3) and a change to this
   document. Being easy to build is not a reason. No dependency or abstraction is added that a
   shipped behaviour does not need.
10. A wrong green is worse than a red: a test or gate that passes while the behaviour is broken is a
    defect of its own.

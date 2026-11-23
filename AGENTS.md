# AGENTS.md

Engineering rules for anyone, human or model, who works on Eulerase. This file is minimal: the
package is being built in steps, and a later step completes it (document roles, setup detail and an
anti-pattern list).

## Setup and gates

```
python -m pip install -e .[dev]
```

Three gates run in this order from the repository root, before every commit:

```
python -m pytest -q -p no:cacheprovider
python -m ruff check .
python tools/check_fresh_export.py
```

The third runs the tests on a copy of the tracked working-tree files (stage new files first). The
core imports neither the window toolkit nor the plotting library, and nothing imports pandas,
scikit-learn or scipy; a test and the linter enforce both.

## Commits

`type(scope): summary`, one concern per commit. The body is for a developer reading `git log`: the
problem with a concrete case, the change, why this design, and which tests pin it. Stage files by
name. Never amend, rebase, force-push or skip the hooks.

## Tests

- A test names the behaviour it guards and fails if that behaviour breaks. It checks an intended
  result against an independent answer, never the code against itself.
- One test per behaviour. Near-duplicates become one parametrized test. No test pins wording or a
  constant.
- No test needs a display, the network, or the user's home directory.
- A wrong green is worse than a red.

## Picking up the work

- **One task, one commit**, or the commits its steps name. Never gather several tasks into one
  commit: a change that cannot be rolled back or reviewed alone is a defect of process.
- **Lean tests.** Write the tests a task names. Add one only where an intended behaviour would
  otherwise have no test, and say so. More tests are not more safety: a test that pins wording,
  repeats another, or checks the code against itself is removed.
- **An unclear instruction is not a conflict.** Resolve it from the section of the design the task
  cites, say in the task's report which reading you took and why, and carry on. Stop and ask only
  for a real contradiction between two rules, or where the choice would change a number, a label or
  a file format a user sees.
- **Carry, then repair.** The verified code is moved as it is; each defect is repaired in its own
  step with its failing test first. Do not rewrite what works.
- **Say what was run.** A report lists the commands that were run with their output, and separately
  what was not run.

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
   evidence of anything in an organism.
9. A parked feature returns only by passing its truth tests and a change to the design. Being easy
   to build is not a reason. No dependency or abstraction is added that a shipped behaviour does not
   need.
10. A wrong green is worse than a red: a test or gate that passes while the behaviour is broken is a
    defect of its own.

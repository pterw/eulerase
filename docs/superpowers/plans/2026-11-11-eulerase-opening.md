# Eulerase v0.1: opening the repository Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Write what a reader and a contributor need on the first day (README, CI workflow, pull request and issue templates, CONTRIBUTING, AGENTS.md), then open the repository, private, with the core, the command line and the window in it, and watch its first CI run.

**Architecture:** Documents and one workflow file. The workflow, the templates and CONTRIBUTING touch nothing the core or the window touches, so they are written beside the other plans' waves. README and AGENTS.md are shared files that the first plan's task A creates in a minimal form; they are completed here, one writer at a time, after the window plan's last task has landed, from facts read out of the landed code. The opening itself is done by the controller, never by an implementer, because it pushes.

**Tech Stack:** GitHub Actions, Markdown, `gh`, git; pytest for the one test of this plan.

**Spec:** `docs/superpowers/specs/2026-11-09-eulerase-v0.1-design.md`, sections 1, 1.1, 1.2, 5, 6.4, 6.6, 7, 8.1, 9, 11, 12 and 12.1 (the authority; where this plan and the spec disagree, the spec wins). Section 9, "Opening the repository", is what this plan carries out. The order of the three plans is: the first plan (`docs/superpowers/plans/2026-11-09-eulerase-v0.1-core.md`, the core and the command line, task ids `A` to `AD`), then the window plan (`docs/superpowers/plans/2026-11-11-eulerase-v0.1-window.md`, the presenter and the window, ids `WA` to `WO`), then this plan's serial tasks and the opening. The ids of the other two plans are used below only in the Dispatch table.

Sections: Global Constraints=brief, Dispatch=reference, Frozen interfaces=brief, File structure=reference, Opening (controller)=reference, What stays open=reference, What cannot be checked here=reference

## Global Constraints

- **The spec is the authority.** A task that cannot be done as the spec says returns NEEDS_CONTEXT; it does not improvise. Parked and dropped features (spec section 3) are not built, described as built, or promised.
- **Only what is shipped.** The opening is the core, the command line and the window. The Windows executable, its measurements and the release report stay out (see "What stays open"). Every document says only what exists at the head it is written for. A document describes the window from the landed code, never from this plan or from memory; no task of this plan opens a window. The README shows a screenshot only if the window plan's hands-on pass produced one (`docs/assets/window.png` is tracked); if that pass was skipped, the README has no image and says in one plain sentence that the window's geometry has not had its hands-on pass. No document says "tested on Linux", "verified" or "supported on Linux" as a result: the workflow is described, a green run is reported by the controller.
- **Independence and the boundary.** Eulerase reads no genotype files, looks up no variant identifiers and interprets nothing for a person (spec 1.1). Documents name no other project except where spec 1.1 does so for the boundary; they do not describe the original explorer or where it lives.
- **Public repository.** No tracked file, commit message or test holds a user-folder path, a drive-letter path, a personal name other than the licence holder, an email address, or a path into any workspace. The source being carried is named only as `$EULERASE_SOURCE`, and only in a plan. This plan is itself a tracked file. The repository opens private and is made public later, once it is packaged; every rule here is kept as if it were public already. Gate lines are relative to the repository root.
- **Three claims kept apart** (spec 1.2): a trail is correct for the graph; a contig is a valid path in the assembly graph built from the k-mers given; a sequence exists in an organism. Eulerase makes the first two and never the third. A contig is not a validated assembly.
- **Multiplicity is not depth** (spec 5.1). Documents say "multiplicity" or "times seen", never "coverage" or "depth", except for simulated reads where the term is correct. The filters (minimum multiplicity, tip clipping) act as error filters only when the input is a set of reads; on a single assembled sequence a low count marks a rarer k-mer, not an error; tip clipping can remove the true end of a sequence. The only two simplifications any document names are those two.
- **Facts come from the landed code.** A document that quotes a help text, a limit, an example or a rule reads it from the tree at its BASE and runs it. It does not copy it from this plan or the spec. Where the landed code and the spec disagree about a fact a document states, the document states what the code does, the report records the difference with both quotations, and the task carries on. It returns NEEDS_CONTEXT only if the difference is a number, a label or a file format that a user sees.
- **Who runs the gates.** A `Mode: code` task of this plan does not run the three gates: it runs only what its change touches (for RA to RE their own checks, since they add no test file; for RH its own test file). The gate-runner runs the gates at landing. The `Mode: serial` tasks (RF, RG) run the commands their documents quote, and the gate-runner follows.
- **No finding, task or review ids in documents** (README, CONTRIBUTING, AGENTS.md, templates). Anti-patterns are written in words.
- **No Dependabot and no scheduled run at opening.** Decided: dependencies are pinned exactly, and both are added later only if they earn their keep (see "What stays open"). There is no CD: nothing is published in v0.1, so there is no deployment, no release job and no secret.
- **Checks that are not tests.** A document is checked by its reviewer against the spec and the landed code, not by a test that pins its wording. The one test of this plan (task RH) compares two documents' commands and fails when they drift.
- **Commits.** `type(scope): summary` as each task names. One concern each. The body explains, for a developer reading `git log`, what the file is for and why it says what it says. No finding, task or review ids, no brief or untracked file cited, no attribution trailer. Stage files by name; never `git add -A`; never amend, rebase, reset, stash, push or `--no-verify`. Wrap at 72.
- **Shell.** Never `sed -i`. A denied tool call is reported, never retried another way. Bound every command with a timeout. Scratch files go in a `mktemp -d` directory outside the repository.
- **Shared files** (`README.md`, `AGENTS.md`, `pyproject.toml`, `LICENSE`, `.gitignore`): in this plan only tasks RF (README) and RG (AGENTS.md) edit them, both `Mode: serial`.

## Dispatch

Frozen: the section "Frozen interfaces" below is part of every task brief. An implementer who needs a name or shape changed returns NEEDS_CONTEXT.

Order of the three plans (owner, 2026-11-11: "the tool" is "the first iteration of the functioning Debruijn graph with the tktinker GUI"): the first plan, the window plan, then this plan's serial tasks and the opening. The first plan's "Left for the second plan" and the window plan's "From the first plan's list" give the README and AGENTS.md in full to this plan, window parts included. This plan takes all of it except the Windows executable, its measurements and the release report, which stay out.

"Starts when" names tasks of the other plans, by their task ids. The tasks of this plan that are `Mode: code` need only the first plan's task A landed (the `.github` files and CONTRIBUTING.md touch nothing the core or the window touches), so they run beside everything. The two `Mode: serial` tasks need the window plan's last landing, WO, and the controller's hands-on pass of that plan run or recorded skipped, so the README knows whether a screenshot exists.

| Wave | Tasks | Starts when | Notes |
|---|---|---|---|
| 1 | RA, RB | A landed | `.github/workflows/gates.yml`, pull request template |
| 2 | RC, RD | A landed | issue templates; wave 2 only caps this plan at two code tasks at once |
| 3 | RE | A landed | `CONTRIBUTING.md`; after RA to RD, whose files it points at |
| 4 | RF | WO landed, the window plan's hands-on pass run or recorded skipped, RE landed | serial; `README.md` |
| 5 | RG | RF landed | serial; `AGENTS.md` |
| 6 | RH | RG landed | the workflow-matches-documented-gates test |

Agent budget. The ceiling for all plans running together is 12 agents, counted as 3 per code task in flight (the task, its nested gate-runner, its reviewer). This plan's code tasks run at most two at a time, so it adds at most 6 agents (waves 1 and 2); wave 3 adds 3, waves 4 and 5 add 1 (serial tasks), wave 6 adds 3. The controller dispatches a task of this plan only when the count of every plan running stays at 12 or below. It never starts one while the first plan's waves 6, 7 or 8 (four code tasks, 12 agents) run, and while the first plan's wave 3 (three tasks, 9 agents) or one of the window plan's waves 2 and 3 (three tasks, 9 agents) runs it starts at most one; while a window wave of two tasks (6 agents) runs it starts at most two. A task that cannot start waits for the next wave boundary of the other plans. The first plan's and the window plan's tasks never wait for this plan.

Landing order: by wave, in id order within a wave (RA, RB, RC, RD, RE, RF, RG, RH). Landings are one at a time on `main`, interleaved with the other plans' landings at whatever point they happen: the controller applies the patch, stages the task's files by name, runs the three gates in the order of the Gates block of Frozen interfaces, commits with the task's subject and dispatches a gate-runner at that HEAD before the next landing. The tasks of this plan touch no file the other plans' tasks touch, except that RF and RG edit the README.md and AGENTS.md that the first plan's task A created, after both other plans have landed; so no conflict is expected, and any conflict returns NEEDS_CONTEXT.

The Gates block (the same three gates as the first plan; the third is added when the first plan's task A lands) is in "Frozen interfaces" below, so every task brief holds it. It is the controller's gate list for the gate-runner.

After RH lands and every task of all three plans is landed, the controller runs "Opening (controller)". No implementer pushes.

## Frozen interfaces

**Files and subjects.**

| Task | File | Commit subject |
|---|---|---|
| RA | `.github/workflows/gates.yml` | `chore(ci): run the gates on Windows and Linux; docs-only PRs skip them` |
| RB | `.github/PULL_REQUEST_TEMPLATE.md` | `docs(github): add a pull request template` |
| RC | `.github/ISSUE_TEMPLATE/bug_report.md` | `docs(github): add a bug report template` |
| RD | `.github/ISSUE_TEMPLATE/feature_request.md` | `docs(github): add a feature request template` |
| RE | `CONTRIBUTING.md` | `docs: add CONTRIBUTING` |
| RF | `README.md` | `docs(readme): describe what is shipped` |
| RG | `AGENTS.md` | `docs(agents): complete the engineering rules` |
| RH | `tests/test_ci_gates.py` | `test(ci): check the workflow runs the documented gates` |

**The Gates block** (read by tasks RB, RG and RH, and by the controller for the gate-runner), exactly:

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
```

The third gate reads tracked files from the working tree: a new file must be staged by name before it runs, or the export lacks it. A document-only task still passes through all three at landing.

**The workflow (task RA).** One file, workflow name `Gates`. Trigger, exactly:

```yaml
on:
  pull_request:
    paths-ignore:
      - docs/**/*.md
      - CONTRIBUTING.md
      - .github/ISSUE_TEMPLATE/*.md
      - .github/PULL_REQUEST_TEMPLATE.md
  push:
    branches: [main]
```

These lines are the **frozen trigger block**: everything from the `on:` line to the line before `permissions:` (which comes next, then an optional `concurrency:`, then `jobs:`), with no comment, blank line or other key inside it. Task RH's test requires the block in the workflow file to equal these lines exactly, after stripping trailing whitespace and line endings; an added, removed or changed trigger line is a flaw.

Documentation-only pull requests skip the gates, and nothing else does. The `paths-ignore` list holds only the four entries above, one per line, unquoted: Markdown files that packaging and the tests never read. They name `.md` files, not whole directories, because a directory pattern would also skip a Python file added under `docs/`, and `python -m ruff check .` lints it; a pull request that adds an image such as `docs/assets/window.png` or any other non-Markdown file under `docs/` therefore runs the gates, which is the safe direction. `README.md` and `LICENSE` are not on it, because a build can read both (a `readme` or licence entry in `pyproject.toml`), so a change to either runs the gates. `AGENTS.md` is not on it (the workflow test reads its gates block), and neither are `src/`, `tests/`, `tools/`, `pyproject.toml` or `.github/workflows/`. GitHub evaluates path filters against a bounded number of changed files, and a skipped pull request is still covered when it reaches `main`, because every push to `main` runs every gate. A pull request that touches any path not on the list runs every gate on both systems, because a workflow skips only when every changed path matches the list; the one exception is a pull request that changes more files than GitHub's filter inspects, which only the push to `main` then covers (a residual risk, named in the final report). `push` to `main` has no path filter: every push to `main` runs every gate. There is no per-module test selection: the suite is small, and a selector that guesses which tests a change reaches is where a false green comes from.

`permissions: contents: read` at the top. One job with id `gates`, `name: gates (${{ matrix.os }})`, `runs-on: ${{ matrix.os }}`, `timeout-minutes: 20`, `strategy.fail-fast: false`, matrix `os: [windows-latest, ubuntu-latest]`. No `defaults:` key and no `shell:` key: each step uses the runner's default shell (PowerShell on Windows, bash on Linux), because the install line and the three gate lines run unchanged under both (the extras are quoted, so neither shell treats `[dev]` as a pattern). Because no `run:` key has an empty value, every `run:` line of the file is a step's command. Steps in order: `actions/checkout` with `persist-credentials: false`; `actions/setup-python` with `python-version: "3.13"`, `cache: pip` and `cache-dependency-path: pyproject.toml` (the repository has no requirements file, so both inputs are needed); then the install and the three gates, each a step with a **single-line** `run:` (no `|` or `>` block, so the lines can be compared). The allowed keys are fixed: a step that has a `run:` has only `name` and `run`, so none has `if`, `continue-on-error`, `working-directory`, `shell`, `env`, `timeout-minutes`, `uses` or `with`; the job has only `name`, `runs-on` (equal to `${{ matrix.os }}`), `strategy` (with `fail-fast` and a `matrix` that holds only `os`, the set `{windows-latest, ubuntu-latest}`, with no `include` or `exclude`), `timeout-minutes` and `steps`, so it has no `continue-on-error` or `if`. The `run:` lines, byte for byte and in this order (the first is the install, not a gate):

```
python -m pip install -e ".[dev]"
python -m pytest -q -p no:cacheprovider
python -m ruff check .
python tools/check_fresh_export.py
```

A `concurrency` group on workflow and ref with `cancel-in-progress: true` is allowed. No `env:` value that hides a difference between the two systems (no forced encoding, no forced backend): a red on one system is a finding. No deployment, release or secret. Actions are pinned to a full commit SHA with the version in a comment if the implementer can verify the SHA offline from a file in the repository; otherwise to the major tag (`actions/checkout@v4`, `actions/setup-python@v5`), and the report says which and why. A one-line comment in the file says there is no CD because nothing is published in v0.1.

**AGENTS.md gates block (task RG; read by task RH).** Under a heading `## Gates`, one fenced block whose info string is `gates`, with the same three items, in the same order, and the same keys as the Gates block above (`name`, `run`, `pass`, `quote`), copied from it with `diff`. Each `run:` value is on one line. This block is the repository's statement of the gate commands; the workflow must equal it.

**Pull request template (task RB).** Headings, in order: `## Summary` (what changed and why), `## How it was checked` (the head SHA, then one line per gate in the order above with its last output line quoted: pytest, ruff, fresh-export; and the CI run link once it is green, or the words "CI skipped: only documentation paths changed" when that is so), `## Not run` (anything not run or not checked, said plainly; "nothing" is allowed only if true). Comments in the file say: quote output, do not paraphrase it; a gate not run is listed under "Not run". No attribution line is asked for.

**Issue templates (tasks RC, RD).** Markdown templates with front matter (`name`, `about`, `title`). `bug_report.md` asks for: the exact command or the input given, the expected result, the actual result with the output copied, the platform and Python version, and `eulerase.__version__` (read with `python -c "import eulerase; print(eulerase.__version__)"`); it warns that the issue is public and that sequences pasted into it are public. `feature_request.md` asks for the problem before the solution; before filing, the reporter checks spec section 3 (the table of kept, parked and dropped features) and section 8.3 (the truth tests a parked feature must pass to return): a parked feature is not requested as if new, and a request for one names the truth test it would pass. It states the boundary of spec 1.1 (no genotype files, no variant lookup, no interpretation for a person). Links to the spec use its path in the repository.

**CONTRIBUTING.md (task RE)** has these sections: what the project is in one sentence; setup (Python 3.13, a virtual environment, `python -m pip install -e ".[dev]"`); the gates (the three commands, their order, that the third reads tracked files so new files are staged first); the commit format (`type(scope): summary`, a body that explains the problem, the change, why this design and which tests pin it, one concern per commit, no ids, no attribution trailer); the test rules of spec 8.1 in short; how a change to behaviour starts (the spec is the authority, a departure from it is made by changing the spec in a commit of its own with the reason, a parked feature returns only through its truth tests); the pull request (the template, the gates quoted, what was not run); a line that CI runs the three gates on Windows and Linux for every push to `main` and for every pull request that changes anything but documentation (the four paths of the workflow's `paths-ignore`, which hold no code, test, tool, `pyproject.toml`, `README.md`, `LICENSE`, `AGENTS.md` or workflow), so a documentation-only pull request shows no CI run, and that nothing is published so there is no release process yet; a pointer to AGENTS.md for the rules an agent or model follows. The test rules include that no test opens a window.

**README.md (task RF)** has: what Eulerase is, in the spec's words for section 1, with the three claims kept apart; the boundary (spec 1.1); install from a clone (nothing is published) and the gates, with the install line quoted as `python -m pip install -e ".[dev]"` (the same in README, CONTRIBUTING.md and the workflow, so a reader who copies it meets no glob in bash); the window: how to start it (`eulerase` with no arguments; `--gui` with input), the layout of spec section 7 as the code builds it (input and settings, the graph with the step controls under it, the five tabs Edges, Nodes, Trails, Contigs and Assembly, the status bar, exports in the File menu), the examples it offers, the Stop control, and the one line the command prints when tkinter is missing; the screenshot `docs/assets/window.png` if that file is tracked, otherwise no image and one plain sentence that the window's geometry has not had its hands-on pass; the command line, with each option group, the exit statuses (0, 1, 2) and the examples; the limits table (spec 6.4: three limits, the k-mer limit, the drawing limit and the table-row limit, each with what happens above it; the minimum ORF length and the default k are not limits and are stated where the ORF rule and the default k are given); what it does not do; what a multiplicity is and the A3 statements about the filters; the statement of the ORF rule (the sentence `orfs.orf_rule_text()` returns) and what the `orfs` column counts; the sentence that the GenBank writer sets its date field to a fixed value, `01-JAN-1980`, so files are byte-identical across runs; the determinism promise in one short paragraph (spec 6.3: text outputs identical; figure data reproducible; image bytes not promised); the licence; a "Status" section saying that the Windows executable is not built, so the tool is run from a clone, and that the repository is private until the tool is packaged.

**Facts RF reads from the landed code** (each quoted in the report with the command that produced it): the output of `eulerase --help`, which includes the `--gui` text; the names of the window's examples, the tab titles, the menu entries with their shortcuts, the Stop control and the confirmation threshold, read from the source of `eulerase.gui.presenter` and `eulerase.gui.app` (read, never run: no task opens a window); the sentence the command prints when tkinter is missing, from `tests/test_cli_window.py`; whether `docs/assets/window.png` is tracked (`git ls-files docs/assets/window.png`); the values of `eulerase.kmers.MAX_KMERS`, `eulerase.render.MAX_DRAW_NODES`, `eulerase.analysis.MAX_TABLE_ROWS`, `eulerase.orfs.MIN_ORF_AA`, `eulerase.analysis.DEFAULT_K`; `eulerase.orfs.orf_rule_text()`; the Contigs table columns from `eulerase.analysis.columns("Contigs")` and how `orfs` is computed; the date literal in `eulerase.io`; and every `main([...])` invocation in `tests/test_cli*.py` (the only examples the README may use; the `--gui` ones in `tests/test_cli_window.py` run against a fake window, so the README may show them as commands, not as output).

**AGENTS.md (task RG)** has these sections: document roles (spec, plans, README, CONTRIBUTING, AGENTS.md: what each is for and which wins); setup (with the install line quoted as `python -m pip install -e ".[dev]"`); `## Gates` (block above); commit rules; test rules (spec 8.1); the window (the presenter, `eulerase.gui.presenter`, holds every decision and imports neither tkinter nor matplotlib; `eulerase.gui.app` is the view and nothing else; no test opens a window or constructs `tk.Tk()`; a test that needs tkinter skips with a stated reason; the hands-on pass is a person's, not a test; the palette is defined once, in `eulerase.render`; only the window's own thread touches a widget; each statement read from spec sections 5 and 7 and checked against the source); the anti-pattern list; the rules of spec 12.1; the ten invariants of spec section 9, copied in the spec's words and numbered as there. The anti-pattern list is seeded from spec section 4 and written in words, without finding ids. It holds at least: a default that removes data (an automatic cutoff); a filter that compares nothing and depends on input order; ordering taken from hashing or from a set; a label or option name that says "coverage" or "depth" for a count of occurrences; a check that reads the code against itself; a test that pins wording or a constant; output taken from the clock, the host or a path; text files written with the platform's line endings; a format guessed from a file name; a path that does not exist read as a sequence; an option accepted and ignored; a failure that ends in a traceback after some output was written; an importing side effect; a table whose header and value come from two definitions; a failed step that leaves an old result visible; a count reported as depth, copy number or evidence of anything in an organism; a parked feature built because it was easy; a decision made in the window instead of the presenter; a test that opens a window; a colour literal outside the palette constants; a widget touched from a worker thread. It does not name the original explorer's location.

## File structure

| File | Responsibility | Task |
|---|---|---|
| `.github/workflows/gates.yml` | the three gates on Windows and Linux | RA |
| `.github/PULL_REQUEST_TEMPLATE.md` | what a pull request states | RB |
| `.github/ISSUE_TEMPLATE/bug_report.md`, `feature_request.md` | what an issue states | RC, RD |
| `CONTRIBUTING.md` | how to work on the project | RE |
| `README.md` | what the tool is and does, for a reader | RF |
| `AGENTS.md` | the engineering rules for any agent or model | RG |
| `tests/test_ci_gates.py` | the workflow equals the documented gates | RH |

### Task RA: The CI workflow

Wave: 1 · Touches: .github/workflows/gates.yml · After: none · Mode: code

Behaviour: every push to `main`, and every pull request that changes anything but documentation (the four paths of Frozen interfaces), runs the three gates, in the first plan's order, on `windows-latest` and `ubuntu-latest`, with Python 3.13. Spec 9 (Gates), 12 step 1. Shape only (never copy names or product text) from the author's other repository: a `Gates` workflow with `pull_request` and `push` to `main`, `permissions: contents: read`, a concurrency group, a job timeout, `persist-credentials: false`.

- [ ] **Step 1: Write the workflow** exactly as "The workflow (task RA)" of Frozen interfaces says: the frozen trigger block copied line for line, permissions, one matrix job, the four single-line `run:` lines in order, no `defaults:` or `shell:` key, and on each step with a `run:` only the keys `name` and `run`. Comments (all outside the trigger block, which holds none): one line that a green run means the same as a green local run, one line that there is no CD because nothing is published in v0.1, one line that only documentation paths are skipped on pull requests and `main` is never filtered.
- [ ] **Step 2: Pins.** Look for a full commit SHA of each action in files of the repository that already hold one (grep `.github/` and `docs/`). If you cannot verify a SHA offline, use `actions/checkout@v4` and `actions/setup-python@v5` and say in the report that the major-tag pin is a known weakness to be raised when the owner chooses to add automated updates. Do not invent a SHA.
- [ ] **Step 3: Check offline.** Quote: the four `run:` lines (`grep -n "run:" .github/workflows/gates.yml` prints exactly four lines, each with a value on the same line); that no `run:` value is a block (`grep -n "run: [|>]"` prints nothing); that the file has no `defaults:`, `shell:`, `working-directory:`, `continue-on-error` or `if:` key (`grep -n -E "defaults:|shell:|working-directory:|continue-on-error|if:"` prints nothing); that the two setup-python cache inputs are present (`grep -n -E "cache: pip|cache-dependency-path: pyproject.toml"` prints two lines); that the trigger block (`grep -n -B1 -A12 "^on:"` quoted) equals the frozen trigger block, with exactly the four `paths-ignore` entries (`docs/**/*.md`, `CONTRIBUTING.md`, `.github/ISSUE_TEMPLATE/*.md`, `.github/PULL_REQUEST_TEMPLATE.md`; the broad `docs/**` is a defect) and `push` as `branches: [main]` alone; that the file parses as YAML if a parser is installed (`python -c "import yaml"` may be absent; if so say "not parsed" and do not install one). Say that no run on a hosted runner was made.
- [ ] **Step 4: Commit** with the subject in Frozen interfaces. Stage the file by name first. This is a `Mode: code` task: it does not run the three gates (it adds no test file; its own checks are Steps 1 to 3); the gate-runner runs them at landing.

Sweep terms: none.

### Task RB: Pull request template

Wave: 1 · Touches: .github/PULL_REQUEST_TEMPLATE.md · After: none · Mode: code

Behaviour: a pull request states what changed and why, how it was checked with the gate output quoted, and what was not run (spec 8.1 "say what was run"; spec 12.1). Shape only from the author's other repository: a summary, a validation list with one line per gate, comments that tell the author what to quote.

- [ ] **Step 1: Write the template** as "Pull request template (task RB)" of Frozen interfaces says. It names no task, finding or plan.
- [ ] **Step 2: Check** that the three gate lines carry exactly the commands of the Gates block (copy them from the Gates block of Frozen interfaces) and that the file holds no path outside the repository. Quote the file's headings (`grep -n "^## "`).
- [ ] **Step 3: Commit** with the subject in Frozen interfaces. The three gates are not run here (no test file is added; Step 2 is the task's check); the gate-runner runs them at landing.

Sweep terms: none.

### Task RC: Bug report template

Wave: 2 · Touches: .github/ISSUE_TEMPLATE/bug_report.md · After: none · Mode: code

Behaviour: a bug report holds what a maintainer needs to reproduce it with no further questions, and warns that the issue is public.

- [ ] **Step 1: Write the template** as "Issue templates (tasks RC, RD)" of Frozen interfaces says (front matter `name: Bug report`, `about`, `title: "bug: "`, no labels since none exist at opening).
- [ ] **Step 2: Check** that the version command it gives runs and prints the version (run it, quote the line), and that front matter parses as `---` delimited lines with the three keys.
- [ ] **Step 3: Commit** with the subject in Frozen interfaces. The three gates are not run here (no test file is added; Step 2 is the task's check); the gate-runner runs them at landing.

Sweep terms: none.

### Task RD: Feature request template

Wave: 2 · Touches: .github/ISSUE_TEMPLATE/feature_request.md · After: none · Mode: code

Behaviour: a feature request starts from a problem, sends the reporter to the table of parked and dropped features and to the truth tests, so a parked feature is not requested as if new, and states the boundary. Spec 3, 8.3, 1.1 and invariant 9.

- [ ] **Step 1: Write the template** as "Issue templates (tasks RC, RD)" of Frozen interfaces says (front matter `name: Feature request`, `about`, `title: "feature: "`). Read spec section 3 and section 8.3 and name the parked features by their spec names (the automatic error model, contig binning, comparison of contigs with a reference), the dropped ones by section, and the truth tests by section number; do not restate the tests.
- [ ] **Step 2: Check** that every spec path and section number it cites exists (`grep -n "^## 3\.\|^### 8\.3" docs/superpowers/specs/*.md`), quoting the output.
- [ ] **Step 3: Commit** with the subject in Frozen interfaces. The three gates are not run here (no test file is added; Step 2 is the task's check); the gate-runner runs them at landing.

Sweep terms: none.

### Task RE: CONTRIBUTING

Wave: 3 · Touches: CONTRIBUTING.md · After: RA, RB, RC, RD · Mode: code

Behaviour: a contributor can set up, run the gates, write a commit and open a pull request from this file alone, and knows that a change to behaviour starts with the spec. Spec 9 (CONTRIBUTING), 8.1, 12.1.

- [ ] **Step 1: Read** the repository's `AGENTS.md` (the minimal one from the first plan's task A), `pyproject.toml` and `tools/check_fresh_export.py` for the commands and the dev dependency group. The commands in the file are the ones the workflow and AGENTS.md carry.
- [ ] **Step 2: Write** `CONTRIBUTING.md` as "CONTRIBUTING.md (task RE)" of Frozen interfaces says. It points at the pull request template and the two issue templates by path. It says nothing about the executable or a release, and about the window only that no test opens one (AGENTS.md has the rest). Its install line is the quoted form of Frozen interfaces.
- [ ] **Step 3: Check.** Run the setup commands it gives in a scratch virtual environment outside the repository if one can be made without network access; if not, say "setup not run" and do not claim it works. Check that the three gate commands it gives equal the Gates block of Frozen interfaces (`diff`, quote the empty output). This is a `Mode: code` task: it does not run the gates; the gate-runner runs them at landing.
- [ ] **Step 4: Commit** with the subject in Frozen interfaces.

Sweep terms: none.

### Task RF: README in full for what is shipped

Wave: 4 · Touches: README.md · After: RE · Mode: serial

Behaviour: a reader learns what Eulerase is, what it does not do, how to install and run it, and sees examples that were run. Every statement is true of the tree at BASE. Spec 9 (README), 1, 1.1, 1.2, 5.1, 6.4, 6.6, 7, B7, A3. Replaces the minimal README of the first plan's task A.

- [ ] **Step 1: Read the facts** listed under "Facts RF reads from the landed code" in Frozen interfaces, by running the commands named there (the window facts by reading the source, never by starting the window), and keep the output for the report. Read the current `README.md` and keep its gate lines. Its install line is the first plan's unquoted `python -m pip install -e .[dev]`; the new README replaces it with the quoted form `python -m pip install -e ".[dev]"`, and Step 6 runs the quoted form.
- [ ] **Step 2: Choose the examples.** List every `main([...])` call in `tests/test_cli*.py` that a reader would recognise as a use (a typed sequence, a set of k-mers, a file, `--metrics`, `--both-strands`, `--csv`, `--gfa`, a failure with exit 1 or 2). The README uses only these, one per kind that exists. A kind with no test is left out and named in the report; do not add a test in this task.
- [ ] **Step 3: Write `README.md`** as "README.md (task RF)" of Frozen interfaces says. The limits table has the three limits read in Step 1 (`MAX_KMERS`, `MAX_DRAW_NODES`, `MAX_TABLE_ROWS`), each with what happens above it (spec 6.4: refused before any work with the limit in the message; the figure says why it is empty; a cut table says how many rows of how many; an export always carries every row). `MIN_ORF_AA` is stated beside the ORF rule and `DEFAULT_K` beside the default k, not in the table. If the first k-mer limit differs from 500,000, the README states the value in the code. The A3 statements about `--min-count` and `--clip-tips` use the help text's meaning, not new claims. The window section follows "README.md (task RF)" of Frozen interfaces: it states only what the source and the tests show, with the source line for each sentence (Step 5). If `git ls-files docs/assets/window.png` lists the file, the README shows it with a link in the repository (relative path), and says nothing about how the window looks beyond what the image and the code show. If it lists nothing, the README has no image and no link to one, and carries the single sentence that the window's geometry has not had its hands-on pass. No claim about speed, size, fit on a screen or platforms tested.
- [ ] **Step 4: Run every example command** of the README from the repository root. The window's own commands (`eulerase`, `eulerase --gui ...`) are not run: they open a window, and no task of this plan does that. They are marked "not run" in the report; their tests in `tests/test_cli_window.py` run against a fake window. The README shows the installed command `eulerase`; run it as that if the console script is on the PATH, otherwise as `python -c "from eulerase.cli import main; raise SystemExit(main([...]))"` with the same arguments, and say in the report which form was used (headless, `MPLBACKEND=Agg`; any file it writes goes to a scratch directory, never into the repository), compare its output with the lines the README shows, and quote each command and its output in the report. The README shows output exactly as the command printed it. Run every exit-status claim the same way. Mark any command not run.
- [ ] **Step 5: Check the claims.** For each sentence that states a number, a rule or an option, quote the line of code, help text or test that makes it true, in a table in the report (sentence, source). A sentence without a source is removed.
- [ ] **Step 6: Run the commands the README quotes, then commit.** This is a `Mode: serial` task: run the install line and the three gate commands as the README gives them, in order, and quote the last line of each; then commit with the subject in Frozen interfaces. The gate-runner follows at that HEAD.

Sweep terms: "coverage", "depth", the to-do markers of the leftover sweep, the unquoted install line, the old command line sentence of the leftover sweep, the task ids of the first and window plans; none may appear in the file. A link to an image appears only if `docs/assets/window.png` is tracked.

### Task RG: AGENTS.md in full

Wave: 5 · Touches: AGENTS.md · After: RF · Mode: serial

Behaviour: any agent or model that picks the work up finds the document roles, the setup, the gates, the commit and test rules, the anti-patterns, the rules of spec 12.1 and the ten invariants in the spec's words. The gates block is machine-read by the test of task RH. Spec 9 (AGENTS.md, Invariants), 12.1, 4. Replaces the minimal AGENTS.md of the first plan's task A.

- [ ] **Step 1: Read** the current `AGENTS.md`, spec section 9 (the invariants), section 12.1 and section 4 (for the anti-patterns). Keep every commit and test rule the minimal file already states; do not weaken one.
- [ ] **Step 2: Write `AGENTS.md`** as "AGENTS.md (task RG)" of Frozen interfaces says. The invariants are copied from spec section 9 word for word and numbered 1 to 10. Where an invariant cites a spec section, the section number stays. The setup line is the quoted install line. The window section states only what the source shows (`git grep` for the imports of `eulerase.gui.presenter` and the absence of `tk.Tk(` in `tests/`, quoted in the report). The `## Gates` block is the fenced `gates` block of Frozen interfaces, copied unchanged.
- [ ] **Step 3: Check the invariants** by a diff of the ten items against the spec text (extract both, compare with `diff`), quoting the empty output; and the gates block against the Gates block of Frozen interfaces the same way.
- [ ] **Step 4: Check the layout facts** it states (source layout, test layout, tools) against `git ls-files`, quoting the listing it used. It states no fact about files that do not exist.
- [ ] **Step 5: Run the commands the document quotes, then commit.** This is a `Mode: serial` task: run the three gate commands of the `## Gates` block, in order, and quote the last line of each; then commit with the subject in Frozen interfaces. The gate-runner follows at that HEAD.

Sweep terms: "coverage" (allowed only inside the anti-pattern that bans it), the task ids of the first and window plans, finding ids, the unquoted install line.

### Task RH: The workflow runs the documented gates

Wave: 6 · Touches: tests/test_ci_gates.py · After: RG · Mode: code

Behaviour: the gate commands CI runs equal the gate commands the repository documents, in the same order, so a change to one that forgets the other is a red. It also guards the skip: the list of paths a pull request may skip the gates for holds only documentation, and `main` is never filtered. It compares documents with each other; it does not pin wording. Spec 9 (Gates, "in CI and before a commit"), invariant 10.

- [ ] **Step 1: Write the tests** in `tests/test_ci_gates.py`, with three helpers in the same file (the third is named below): `documented_gates(agents_text: str) -> list[str]` (the `run:` values of the fenced block whose info string is `gates` in `AGENTS.md`) and `workflow_runs(workflow_text: str) -> list[str]` (the values of every `run:` line of `.github/workflows/gates.yml`; the frozen workflow has no `defaults:` key, so every `run:` line is a step's command; raises `AssertionError` on a `run:` whose value is `|`, `>` or empty, so a block cannot hide a command). Both files are read with `Path.read_text(encoding="utf-8")`, never as bytes, and every line and value is compared after `strip()`, so a CRLF checkout on Windows leaves no `\r` on a value. Tests: `test_workflow_runs_the_documented_gates_in_order` (read both files from the repository root through `Path(__file__).resolve().parents[1]`; the workflow's runs after its first one, the install line, equal the documented list; the install line equals `python -m pip install -e ".[dev]"` exactly (it is the one `run:` the documented list does not cover, so a second command after `;` or `&&` on it would run unseen); the list has at least the three gates; `workflow_flaws(workflow_text) == []`, the third helper, which makes a gate that cannot fail, a matrix that lost a system and a trigger that changed a red). The test **allows; it does not forbid**: it holds the workflow to the shape the plan freezes, because a list of forbidden keys can always be outrun by a break the list does not name. The third helper first cuts every ` #` comment from the workflow's lines (keeping indentation) and drops blank lines, so a name that sits only in a comment is not seen; the trigger check below reads the lines as they are, before that cut. It lists a flaw for each of these:
  - **Trigger.** The lines from `on:` to the line before `permissions:` or `jobs:`, with trailing whitespace and line endings stripped, are not equal to the frozen trigger block of Frozen interfaces, line for line (an added, removed or changed line, a comment or a blank line inside the block is a flaw; this one equality covers `paths`, `paths-ignore` under `push`, `branches`, `branches-ignore`, `types` and a removed trigger).
  - **Gate steps.** A step that has a `run:` key has any key other than `name` and `run` (`if`, `continue-on-error`, `working-directory`, `shell`, `env`, `timeout-minutes`, `uses`, `with` and all others).
  - **Job.** The job `gates` has a key other than `name`, `runs-on`, `strategy`, `timeout-minutes` and `steps` (`continue-on-error` and `if` among them); its `runs-on` is not exactly `${{ matrix.os }}` (a hard-coded system is a flaw); its `strategy` has a key other than `fail-fast` and `matrix`, or its `fail-fast` line is not exactly `fail-fast: false` (the literal, so `true` and `${{ true }}` are flaws); its `matrix` has a key other than `os`, or an `os:` that is not one line `os: [a, b]` whose two entries, compared as a set, are exactly `{windows-latest, ubuntu-latest}` (a substring test is not enough: a system left in a comment, an `exclude:` entry or an `include:` would pass it).
  - **File.** A top-level key other than `name`, `on`, `permissions`, `concurrency` and `jobs`, or a job other than `gates`.

  `test_oracle_rejects_a_drifted_workflow` is **one** parametrized test over the frozen workflow text with one break each; every case is reported as different, flawed or raises, and the real text is the control. The cases are at least: a gate `run:` line altered; `src/**` added to `paths-ignore`; `branches-ignore: [main]` under `push`; the `push` trigger removed; `types: [closed]` under `pull_request`; `ubuntu-latest` removed from the matrix but left in a comment; `continue-on-error: true` on a gate step; `working-directory: docs` on a gate step; `if: false` on the job. a second command appended to the install line after `; `; `fail-fast: true`; `docs/**` in place of `docs/**/*.md` in `paths-ignore`; a seventh step before the gates. Cases the first plan's shape also needs stay: one gate removed, two gates swapped, an extra `run:` line between the gates, a `run: |` block, a hard-coded `runs-on: ubuntu-latest`, an `exclude:` entry.
- [ ] **Step 2: Red.** Show the test fails on drift without touching tracked files: copy `.github/workflows/gates.yml` and `AGENTS.md` into a scratch directory, change one flag in the copy of the workflow (and, in a second copy, add `src/**` to its `paths-ignore`; in a third, put `branches-ignore: [main]` under `push`), run the helpers on the copies (a short `python -c` over the helpers imported from the test file) and quote that the first differs and that the second and third have a flaw. Then the control: on the real files they are equal and `workflow_flaws` is empty.
- [ ] **Step 3: Green.** Run `python -m pytest tests/test_ci_gates.py -q -p no:cacheprovider` and quote the last line. The file must also pass under the third gate (it reads only tracked files).
- [ ] **Step 4: Commit** with the subject in Frozen interfaces. This is a `Mode: code` task: it runs its own test file (Step 3), not the three gates; the gate-runner runs them at landing, and the file must pass under the third gate. The body says the test compares the workflow with the documented gates and that it does not check what the gates do.

Sweep terms: none.

## Opening (controller)

The controller does this section after every task of all three plans has landed. Implementers never push. Each command is run with a timeout; its output goes to the final report. Scratch files go in a `mktemp -d` directory outside the repository.

1. **Preconditions.** Every task of all three plans (the first plan, the window plan and this plan) is landed, reviewed and gated, and the window plan's hands-on pass has been run or recorded skipped. On `main`'s head the three gates are green (run them once more; quote the last line of each). `git status --short` shows no tracked change, and no untracked file is expected to be part of the repository. `git remote -v` prints nothing (there is no remote yet). `gh auth status` succeeds and its output lists the `workflow` scope (a push that adds `.github/workflows/` is refused without it). A gate that is red or a tracked change is not an owner question: the opening does not begin, and the controller goes back to the task that left it (a fix task, then the gates again). A remote that already exists, `gh auth status` failing or the `workflow` scope missing are waits for the owner (see "Where the opening may wait for the owner" at the end of this section).
2. **Privacy check.** The controller writes the owner's names to a scratch file `names.txt` (given name, surname, account name, the user-folder name of the machine), one per line, outside the repository; they are not written into any tracked file or message. First `test -s names.txt` must succeed: an empty or missing file makes every names command below pass without checking anything, so the controller rebuilds it before going on. Then the positive control: `git grep -n -I -i -F -f names.txt -- LICENSE` must print at least one line (the licence holder's name); if it prints nothing, the names file or the command is wrong: the controller rebuilds the names file and runs the control again (not a wait for the owner). Then, from the repository root:
   - Tracked files, a drive-letter path (a letter that follows a non-letter or the start of the line, so the scheme of a URL does not match): `git grep -n -I -E '(^|[^A-Za-z])[A-Za-z]:[/\\]' -- .`
   - A user-folder segment: `git grep -n -I -E '[Uu]sers[/\\]|/home/[a-z]' -- .`
   - An address: `git grep -n -I -E '[[:alnum:]._+-]+@[[:alnum:]-]+[.][[:alpha:]]+' -- .`
   - A workspace path: `git grep -n -I -E '[.]superpowers|ss[w][/\\]|sd[d][/\\]' -- .`
   - The owner's names: `git grep -n -I -i -F -f names.txt -- .`
   - Every commit message through the same patterns: `git log --format=%B | grep -n -E -e '(^|[^A-Za-z])[A-Za-z]:[/\\]' -e '[Uu]sers[/\\]' -e '[[:alnum:]._+-]+@[[:alnum:]-]+[.][[:alpha:]]+' -e '[.]superpowers'` and `git log --format=%B | grep -n -i -F -f names.txt`.
   - The identities: the commit identity of this repository is the owner's global git identity, as in the owner's other repositories, so no identity is printed and none is a stop. The check is that the set has not grown since the plan was written, as a count that prints only a number: `git log --format='%an <%ae>%n%cn <%ce>' | sort -u | wc -l` prints `1` (the count at the commit this plan was written on; measured, not assumed). A larger number means a task committed under another identity: the controller stops and reports the count, never the identities. This is one of the places where the opening waits for the owner.

   **Pass** means: each command prints nothing, except (a) the licence holder's name in `LICENSE` and the first plan's line that quotes the licence holder line, for the names command only, and (b) hits that are the pattern text itself inside a plan, which the controller reads and lists one by one (the patterns above are written so they do not match themselves; a hit is a defect until read). The identity count passes only at the number above. A hit in a tracked file is fixed by a task before going on, never by an edit made during the opening. A hit in an existing commit message cannot be fixed by a task, and rewriting history is outside this plan: the controller stops, reports which pattern matched and in which commit (its short SHA, never the matched text), and pushes nothing. This is one of the places where the opening waits for the owner.
3. **Sweep for leftovers.**
   - `git grep -n -I -E 'TO[D]O|FIXM[E]|XX[X]' -- .` prints nothing.
   - `git grep -n -I 'indow is not buil[t]' -- . ':!docs/superpowers'` prints nothing: the window is built, so the old command line sentence appears in no README, document, source file or test. Only the plans and the spec, which are point-in-time documents, may quote it; a hit anywhere else is a document that denies a thing the code decides, or a task that left the old line behind; fix it by a task.
   - The install line is quoted everywhere a reader copies it: `git grep -n -I -F 'install -e .[dev]' -- . ':!docs/superpowers'` prints nothing (the quoted form `-e ".[dev]"` does not match it).
   - The first plan's task ids in tracked files that are not plans: `git grep -n -I -E '\b[Tt]ask [A-Z][A-Z0-9]*\b' -- . ':!docs/superpowers/plans'` prints nothing.
   - Finding ids in the documents: `git grep -n -E '\b(A[1-8]|B[2679]|D([2-9]|1[02])|G1|R[1-3])\b' -- README.md AGENTS.md CONTRIBUTING.md .github` prints nothing (read each hit; a word that happens to match is listed, not edited).
   - `git ls-files` lists no file the repository should not hold (no scratch, cache, export or `.env`).
4. **Create the repository and push.** Owner, 2026-11-11: the goal is "the Eulerase repo initial commit open with the README.md, AGENTS.md, the tool, the CI/CD config"; "It can be kept private." So:
   - `gh repo create eulerase --private --source . --remote origin` (no `--push`).
   - `git rev-parse HEAD` gives `<sha>`; confirm it is the head the gates ran on. Then `git push --atomic origin <sha>:refs/heads/main`.
   - `gh repo view eulerase --json url,visibility` and quote the address and that the visibility is private.
   - The repository stays private. Owner, 2026-11-11: "Private for now, but treat as public since I'd be opening it once it is packaged as an executable." The controller never makes it public; that is done later, once the executable exists (see "What stays open").
5. **Watch the first run.** Before any watcher is dispatched, replace the whole `## CI` section of the controller's constraints file (the one that says "None yet" and forbids a watcher) with this, one `ci` block under a heading that contains `CI`:

   ````
   ## CI

   ```ci
   workflow: Gates
   event: push
   jobs: 2
   blocking:
     - 'gates (*)'
   ```
   ````

   `workflow` is the workflow's `name`; `event` is the event of the run (`push` to `main`; the watcher is given the branch separately); `jobs: 2` is the exact count (`gates (windows-latest)` and `gates (ubuntu-latest)`); the one glob makes both blocking. Nothing is advisory, so there are no other keys. Then dispatch a `ci-watcher` in the background with these fields: `WORKSPACE` (the plan's workspace directory, the one holding that constraints file), `HEAD` (the literal 40-character `<sha>` that was pushed, never `HEAD` or a short SHA, which the watcher refuses), `BRANCH` `main`, `PREFIX` (a log prefix matching `[a-z][a-z0-9-]{0,39}`, for example `open-ci1`; a new one for each loop of step 6), and `REPO` (this repository's root). Wait for its summary and quote each job's result. If it answers `NO_RUN`, wait until `gh run list --branch main --json headSha` shows the head, then dispatch again with `RUN_ID`.
6. **If the first run is red.** Linux has never run this code, so a red is expected to be possible. The failing job's log names the failing test or step. Say apart three classes: (a) a red that is the test's fault (it assumed Windows paths, line endings, a locale, a tool that is not on the runner); (b) a red that is the code's fault (the tool behaves differently on Linux); (c) a red that is neither a test nor the code: the workflow file or the runner setup (YAML, an input name, the pip cache setting, a missing token scope, a missing tool on the image). Then, in this order: a fresh implementer fixes it in a code phase (for (a) the failing test first, made correct for both systems without being weakened; for (b) the failing test first, the test stays and the code is fixed; for (c) there is no test to write first, so the change is to the workflow file, with the failing job's log line quoted as the evidence in its place, and the test of task RH still passes); a reviewer reads it; the controller lands it, runs the three gates locally, and pushes by SHA with `git push --atomic origin <new-sha>:refs/heads/main`; a new watcher follows that SHA. At most three such loops, of any class together. After the third red, or when a red is not reproducible or has a cause outside the repository (a runner image, a service), the controller stops and reports what is red and why, without further changes. A red is never made green by skipping or loosening a test or removing a matrix entry.
7. **Final report.** The repository's address; the head SHA; the CI run and its result per job (green, red with the reason, or not finished); the privacy and leftover checks with their output; and everything "not checked": that a pull request with more changed files than GitHub's path filter inspects was never skipped wrongly, the Windows executable and its measurements, the release report, the hands-on items of the window plan that were skipped or not checked (and whether the README has a screenshot), tests skipped for want of tkinter on a runner, the `KC` tag against the GFA text, anything Linux showed, any behaviour on a Python other than 3.13, and that the actions are pinned to major tags if RA could not verify a SHA.

**Where the opening may wait for the owner, and nowhere else.** (1) `gh auth status` fails, lacks the `workflow` scope, or a remote already exists. (2) `gh repo create` fails because the name is taken or the permission is missing; there is no fallback name. (3) A private string in an existing commit message (step 2). (4) The fix loop's third red, or a red with a cause outside the repository (step 6). (5) The identity count of step 2 above the stated number. Every other case has a default in this plan: a hit in a tracked file is fixed by a task, a red gate or a tracked change sends the controller back to a fix task, a failed positive control or an empty names file is fixed by rebuilding the names file, a missing network for the setup check is reported as "setup not run", and an action without a verifiable SHA is pinned to its major tag and reported.

## What stays open

- **The Windows executable, its measurements and the release report** (spec 12 step 8, spec 9): a later release step; nothing here builds them. The repository opens private and is made public later, once it is packaged. Owner, 2026-11-11: "Private for now, but treat as public since I'd be opening it once it is packaged as an executable."
- **The window's hands-on items** that the window plan's pass skipped or did not check, and the screenshot if that pass made none. The README says so in one sentence and is amended when the pass is done.
- **The `KC` check.** The definition of `KC`, and the use of a lower-case `st` tag on the GFA header line, are not confirmed against the text of the GFA 1 specification. It needs the owner's approval for network access. It is a release item (spec 11), not an opening item.
- **Dependabot and a scheduled run.** Decided not to add them at opening (dependencies are pinned exactly). Revisit at the release, when the report records a pinned dependency with a published advisory, or when the owner asks; a scheduled run would also be the first way to notice a change in a runner image. Pinning actions to SHAs goes with the same decision if RA could only use major tags.
- **Branch protection and required checks.** Not set by this plan; the controller reports that they are not set. If required checks are set later, a workflow skipped by `paths-ignore` leaves a required check pending, so the owner decides how documentation-only pull requests are merged then.

## What cannot be checked here

- **A hosted run.** The workflow is written and read offline. Its syntax, the runner images and the behaviour on Linux are first shown by the run in Opening step 5. No result for Linux is claimed before it.
- **The pins of the actions.** A SHA cannot be verified without the network; if RA uses major tags, that is reported, not hidden.
- **Document truth beyond the checks.** The README and AGENTS.md are checked by running their commands and by quoting a source for each claim, and by a reviewer. No test pins their wording, so a later change to the code can make a sentence stale; the reviewer of any task that changes a help text, a limit or a rule reads the README for it.
- **The window itself.** The README describes the window from its source and tests; no task here opens one. Its geometry, wrapping and responsiveness are checked only if the window plan's hands-on pass ran.
- **Setup from a clean machine.** RE runs the setup only if a scratch environment can be made without the network.
- **Platform line endings of the checkout.** The first plan's byte tests write files with `\n`; whether a Windows checkout changes the line endings of tracked files is first seen in the hosted run.

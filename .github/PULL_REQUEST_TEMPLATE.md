## Summary

<!-- What changed and why. Say the problem first, then the change. -->

## How it was checked

<!--
Quote the last line of each gate's output. Do not paraphrase it.
A gate you did not run is not listed here: put it under "Not run".
-->

Head SHA:

- pytest: `python -m pytest -q -p no:cacheprovider`
  Last output line:
- ruff: `python -m ruff check .`
  Last output line:
- fresh-export: `python tools/check_fresh_export.py`
  Last output line:

CI run (once it is green):

<!--
If only documentation paths changed, CI does not run. Replace the line above
with: CI skipped: only documentation paths changed
-->

## Not run

<!--
Anything not run or not checked, said plainly. "Nothing" is allowed only if
it is true.
-->

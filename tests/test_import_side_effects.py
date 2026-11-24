"""Importing the package or any of its modules must not touch the interpreter's global state."""

import json
import pkgutil
import subprocess
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

import pytest
from genomes import child_env

# The streams are replaced first, so the start is the same on every platform and is not what a
# reconfigure would set. Prints a JSON list of what the imports changed (empty when clean).
PROBE = """
import importlib, io, json, os, sys

real_out = sys.stdout


def fresh():
    return io.TextIOWrapper(
        io.BytesIO(), encoding="cp1252", errors="backslashreplace", newline="\\r\\n",
        write_through=True,
    )


streams = {"stdout": fresh(), "stderr": fresh()}
sys.stdout, sys.stderr = streams["stdout"], streams["stderr"]


def stream_state(name):
    stream = getattr(sys, name)
    buffer = streams[name].buffer
    buffer.seek(0)
    buffer.truncate()
    stream.write("\\u00b7\\n")
    return [stream is streams[name], stream.encoding, stream.errors, buffer.getvalue().hex()]


path, environ = list(sys.path), dict(os.environ)
states = {name: stream_state(name) for name in streams}
for module in sys.argv[1:]:
    importlib.import_module(module)
changes = []
if sys.path != path:
    changes.append("sys.path changed")
if dict(os.environ) != environ:
    changes.append("os.environ changed: " + ", ".join(
        sorted(k for k in set(environ) | set(os.environ) if environ.get(k) != os.environ.get(k))))
for name in streams:
    if stream_state(name) != states[name]:
        changes.append("sys.%s changed: %s -> %s" % (name, states[name], stream_state(name)))
real_out.write(json.dumps(changes))
"""


def import_side_effects(
    modules: Sequence[str], env: Mapping[str, str] | None = None, cwd: str | None = None
) -> list[str]:
    """Import `modules` in a child interpreter; list what changed in sys.path, streams, environ."""
    done = subprocess.run(
        [sys.executable, "-c", PROBE, *modules],
        capture_output=True,
        text=True,
        timeout=120,
        env=child_env() if env is None else dict(env),
        cwd=cwd,
    )
    assert done.returncode == 0, done.stderr
    return json.loads(done.stdout)


PACKAGE_DIR = Path(__file__).resolve().parents[1] / "src" / "eulerase"
MODULES = ["eulerase"] + [
    f"eulerase.{info.name}" for info in pkgutil.iter_modules([str(PACKAGE_DIR)])
]


@pytest.mark.parametrize("module", MODULES)
def test_importing_the_package_and_every_module_changes_nothing(module, tmp_path):
    assert {"eulerase", "eulerase.render", "eulerase.io"} <= set(MODULES)
    env = child_env()
    env["MPLCONFIGDIR"] = str(tmp_path)
    env["MPLBACKEND"] = "Agg"
    assert import_side_effects([module], env=env) == []


def _bad_module_reports(tmp_path, body: str) -> list[str]:
    (tmp_path / "bad_module.py").write_text(body, encoding="utf-8")
    return import_side_effects(["bad_module"], cwd=str(tmp_path))


def test_oracle_reports_a_module_that_edits_sys_path(tmp_path):
    reported = _bad_module_reports(tmp_path, "import sys\nsys.path.append('/elsewhere')\n")
    assert len(reported) == 1 and reported[0].startswith("sys.path")


def test_oracle_reports_a_module_that_reconfigures_stdout(tmp_path):
    reported = _bad_module_reports(tmp_path, "import sys\nsys.stdout.reconfigure(newline='\\n')\n")
    assert len(reported) == 1 and reported[0].startswith("sys.stdout")

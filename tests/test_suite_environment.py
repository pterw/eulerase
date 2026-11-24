"""The suite's matplotlib setup: headless, and its font cache is kept out of the home directory."""

import os
import tempfile
from pathlib import Path


def test_matplotlib_runs_headless_with_its_cache_in_a_temporary_directory():
    assert os.environ["MPLBACKEND"] == "Agg"
    config = Path(os.environ["MPLCONFIGDIR"]).resolve()
    assert Path(tempfile.gettempdir()).resolve() in config.parents
    assert config != Path.home().resolve() / ".matplotlib"

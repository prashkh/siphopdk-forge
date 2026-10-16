import pathlib
import re

import siphopdk


def test_version():
    assert isinstance(siphopdk.__version__, str)
    assert siphopdk.__version__ == siphopdk.siphopdk().version
    pyproject = pathlib.Path(__file__).resolve().parent.parent / "pyproject.toml"
    if pyproject.is_file():
        match = re.search('^version = "([^"]*)"', pyproject.read_text(), re.M)
        assert match and match[1] == siphopdk.__version__
    for rules in (siphopdk.DRC_RULES, siphopdk.LVS_RULES):
        assert f"v{siphopdk.__version__}" in rules.read_text().splitlines()[1]

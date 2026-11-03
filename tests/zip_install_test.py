"""The rule files are found when the package is installed as a zip archive."""

import subprocess
import sys
import zipfile
from pathlib import Path

import siphopdk


def test_rule_files_are_resources():
    for rules in (siphopdk.DRC_RULES, siphopdk.LVS_RULES):
        assert rules.is_file()
        with siphopdk.rules_path(rules) as path:
            assert path.is_file()
            assert path.name == rules.name
            assert path.read_text() == rules.read_text()


def test_rules_load_from_a_zip_archive(tmp_path):
    """A copy of the package inside a zip archive, which has no directory on disk."""
    package = Path(siphopdk.__file__).parent
    archive = tmp_path / "pdk.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        for file in package.rglob("*"):
            if file.is_file() and "__pycache__" not in file.parts:
                zf.write(file, Path(package.name) / file.relative_to(package))

    code = (
        "import sys\n"
        f"sys.path.insert(0, {str(archive)!r})\n"
        "import siphopdk as pdk\n"
        "assert '.zip' in pdk.__file__, pdk.__file__\n"
        "for rules in (pdk.DRC_RULES, pdk.LVS_RULES):\n"
        "    with pdk.rules_path(rules) as path:\n"
        "        assert path.is_file() and path.read_text()\n"
    )
    subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)

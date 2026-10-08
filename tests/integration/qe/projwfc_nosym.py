"""Apply the documented MgO projection workaround to stdin, then run QE."""

import re
import subprocess
import sys
from pathlib import Path


def disable_symmetry(text):
    modified, count = re.subn(
        r"(?im)^(\s*lsym\s*=\s*)\.(?:true|false)\.",
        r"\g<1>.false.", text,
    )
    if count != 1:
        raise ValueError("Expected exactly one lsym assignment in projwfc input")
    return modified


def main():
    modified = disable_symmetry(sys.stdin.read())
    Path("proj-used.in").write_text(modified)
    print("[CI projection policy] lsym=.false.; see tests/integration/qe/README.md", flush=True)
    return subprocess.run(["projwfc.x"], input=modified, text=True).returncode


if __name__ == "__main__":
    sys.exit(main())

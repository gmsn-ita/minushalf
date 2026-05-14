"""
Mock VASP script called by minushalf execute via minushalf.yaml.

Instead of running a real VASP calculation, this script copies pre-computed
output files from INPUTS/ to the current working directory.

Usage in minushalf.yaml:
    software: VASP
    vasp:
        command: ['python', '/full/path/to/metalcheck.py']
    correction:
        correction_code: v
"""
import shutil
import sys
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

def copy_mock_output():
    """
    Copy all files and folders from INPUTS/ into the
    current working directory (wherever minushalf is running from).
    """

    # Path minushalf is running into
    dest_dir = Path.cwd()

    # INPUTS directory
    source_dir = Path(__file__).parent / "INPUTS"

    if source_dir.exists():
        for item in source_dir.iterdir():
            dest = dest_dir / item.name
            if item.is_file():
                shutil.copy2(item, dest)
            elif item.is_dir():
                if dest.exists():
                    shutil.rmtree(dest)
                shutil.copytree(item, dest)
    else:
        print(f"WARNING: INPUTS directory not found at {source_dir.resolve()}")
    print(f"\nDONE\n")  

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print(f"\n[SIMULATED VASP SCRIPT RUNNING]")

    # Copy the corresponding mock output files
    copy_mock_output()
    sys.exit(0)
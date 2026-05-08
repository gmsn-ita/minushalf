#!/usr/bin/env python3
"""
mock_vasp.py
------------
Mock VASP script called by minushalf execute via minushalf.yaml.

Instead of running a real VASP calculation, this script copies pre-computed
output files from INPUTS/cut_2.XX to the current working directory.

Each time this script is called, it advances through the CUT_SEQUENCE array,
so successive calls copy from different input folders, simulating multiple
VASP runs with different cutoff values.

Usage in minushalf.yaml:
    software: VASP
    vasp:
        command: ['python', '/full/path/to/sanitycheck.py']
    correction:
        correction_code: v

To reset the sequence between test runs, delete the .mock_vasp_state file.
"""
import shutil
import sys
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

def copy_mock_output():
    """
    Copy all files and folders from INPUTS/cut_2.XX into the
    current working directory (wherever minushalf is running from).
    """

    # Folder containing the pre-computed VASP output folders (cut_2.59, cut_2.72, etc.)
    # Lives next to this script.
    inputs_dir = Path(__file__).parent / "INPUTS"

    # Build the exact destination path minushalf expects
    dest_dir = Path.cwd()

    # Get the last part of the path (folder name)
    folder_name = dest_dir.name

    # Extract the value after 'cut_2.'
    if folder_name.startswith("cut_2."):
        cut_value = folder_name.split("cut_2.")[1]
        source_dir = inputs_dir / f"cut_2.{cut_value}"

        # Verify the source folder actually exists
        if not source_dir.exists():
            print(f"[MOCK VASP] ERROR: Source folder not found: {source_dir.resolve()}")
            print(f"[MOCK VASP] Make sure INPUTS/cut_2.{cut_value} exists next to this script.")
            sys.exit(1)

        items = list(source_dir.iterdir())
        # Copy each item from the source folder to the current working directory
        for item in items:
            dest = dest_dir / item.name

            if item.is_file():
                shutil.copy2(item, dest)

            elif item.is_dir():
                # Remove destination folder first if it already exists
                if dest.exists():
                    shutil.rmtree(dest)
                    print(f"[MOCK VASP] Removed existing folder: {item.name}")
                shutil.copytree(item, dest)

        print(f"\nDONE\n")
    else:
        source_dir = inputs_dir / "cut_2.dummy"
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
            print(f"[MOCK VASP] WARNING: cut_2.dummy not found at {source_dir.resolve()}")
        print(f"\nDONE\n")  

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print(f"\n[SIMULATED VASP SCRIPT RUNNING]")
    # Copy the corresponding mock output files
    copy_mock_output()
    sys.exit(0)
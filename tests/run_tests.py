#!/usr/bin/env python3
"""
run_tests.py
------------
Reads the test registry from mhtests.yml and runs each test sequentially.

For each test it:
1. Resolves the absolute path of the test script
2. Writes that path to the minushalf.yaml file
3. Runs 'minushalf execute' inside the test folder
4. Reports pass or fail

Usage:
    python run_tests.py

To add a new test, just add an entry to mhtests.yml — no changes needed here.
"""

import subprocess
import sys
import yaml  
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

# Directory where run_tests.py lives (the tests/ folder)
TESTS_ROOT = Path(__file__).parent.resolve()

# Path to the test registry
REGISTRY_FILE = TESTS_ROOT / "testconfig.yml"

# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

import os
import shutil

def potcar_setup(base_dir="test_farm/VASP"):
    """
    Copies and organizes POTCAR files for VASP test cases.
    """

    potfiles_dir = os.path.join(base_dir, "POTFILES")

    # Helper functions
    def ensure_dir(path):
        os.makedirs(path, exist_ok=True)

    def copy_file(src, dst):
        shutil.copy2(src, dst)

    def cat_files(file_list, output_file):
        with open(output_file, 'wb') as outfile:
            for fname in file_list:
                with open(fname, 'rb') as infile:
                    outfile.write(infile.read())

    # === 1. Carbon ===
    carbon_input_dir = os.path.join(base_dir, "1_carbon", "INPUTS", "cut_2.dummy")
    carbon_mh_dir = os.path.join(base_dir, "1_carbon", "minushalf_potfiles")

    ensure_dir(carbon_input_dir)
    ensure_dir(carbon_mh_dir)

    potcar_c = os.path.join(potfiles_dir, "POTCAR.c")

    copy_file(potcar_c, os.path.join(carbon_mh_dir, "POTCAR.c"))
    copy_file(potcar_c, os.path.join(carbon_input_dir, "POTCAR"))

    # === 2. Metal ===
    metal_input_dir = os.path.join(base_dir, "2_metal", "INPUTS")
    metal_mh_dir = os.path.join(base_dir, "2_metal", "minushalf_potfiles")

    ensure_dir(metal_input_dir)
    ensure_dir(metal_mh_dir)

    potcar_ag = os.path.join(potfiles_dir, "POTCAR.ag")
    potcar_f = os.path.join(potfiles_dir, "POTCAR.f")

    copy_file(potcar_ag, os.path.join(metal_mh_dir, "POTCAR.ag"))
    copy_file(potcar_f, os.path.join(metal_mh_dir, "POTCAR.f"))

    cat_files(
        [potcar_ag, potcar_f],
        os.path.join(metal_input_dir, "POTCAR")
    )

    # === 3. NND ===
    nnd_input_dir = os.path.join(base_dir, "3_nnd", "INPUTS")
    nnd_mh_dir = os.path.join(base_dir, "3_nnd", "minushalf_potfiles")

    ensure_dir(nnd_input_dir)
    ensure_dir(nnd_mh_dir)

    potcar_sb = os.path.join(potfiles_dir, "POTCAR.sb")
    potcar_s = os.path.join(potfiles_dir, "POTCAR.s")

    copy_file(potcar_ag, os.path.join(nnd_mh_dir, "POTCAR.ag"))
    copy_file(potcar_sb, os.path.join(nnd_mh_dir, "POTCAR.sb"))
    copy_file(potcar_s, os.path.join(nnd_mh_dir, "POTCAR.s"))

    cat_files(
        [potcar_ag, potcar_sb, potcar_s],
        os.path.join(nnd_input_dir, "POTCAR")
    )

    print("POTCAR setup completed successfully.")

def load_registry():
    """Load and return the list of tests from mhtests.yml."""
    if not REGISTRY_FILE.exists():
        print(f"[ERROR] Registry file not found: {REGISTRY_FILE}")
        sys.exit(1)

    with open(REGISTRY_FILE, "r") as f:
        data = yaml.safe_load(f)

    return data.get("tests", [])


def run_test(test: dict):
    """
    Run a single test entry from the registry.

    Args:
        test (dict): A test entry from mhtests.yml with keys:
                     name, description, folder, script
    """
    name        = test["name"]
    description = test["description"]
    folder      = TESTS_ROOT / test["folder"]
    script      = folder / test["script"]

    print(f"\n{'='*60}")
    print(f"  TEST : {name}")
    print(f"  DESC : {description}")
    print(f"  DIR  : {folder}")
    print(f"  SCRIPT: {script}")
    print(f"{'='*60}\n")

    # Verify the test folder exists
    if not folder.exists():
        print(f"[ERROR] Test folder not found: {folder}")
        return False

    # Verify the test script exists
    if not script.exists():
        print(f"[ERROR] Test script not found: {script}")
        return False

    # Write minushalf.yaml with the absolute path of the script injected
    minushalf_yaml = folder / "minushalf.yaml"
    minushalf_yaml.write_text(
        f"software: VASP\n"
        f"vasp:\n"
        f"    command: ['python', '{script.resolve()}']\n"
        f"correction:\n"
        f"    correction_code: v\n"
    )

    # Run minushalf execute inside the test folder
    print(f"[INFO] Running 'minushalf execute' in {folder}\n")
    result = subprocess.run(
        ["minushalf", "execute"],
        cwd=folder,
        stdout=sys.stdout,
        stderr=sys.stderr
    )

    # Report result
    if result.returncode == 0:
        print(f"\n[PASS] ✅ {name} passed.")
        return True
    else:
        print(f"\n[FAIL] ❌ {name} failed with return code {result.returncode}.")
        return False


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print(f"\n{'='*60}")
    print(f"  MHTEST - Minushalf Test Runner")
    print(f"  Tests root : {TESTS_ROOT}")
    print(f"  Registry   : {REGISTRY_FILE}")
    print(f"{'='*60}")

    potcar_setup()
    tests = load_registry()
    print(f"\n[INFO] Found {len(tests)} test(s) to run.\n")

    results = {}
    for test in tests:
        passed = run_test(test)
        results[test["name"]] = passed

    # ── Final summary ────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"  SUMMARY")
    print(f"{'='*60}")
    for name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"  {status}  {name}")

    total  = len(results)
    passed = sum(results.values())
    failed = total - passed
    print(f"\n  Total: {total}  |  Passed: {passed}  |  Failed: {failed}")
    print(f"{'='*60}\n")

    # Exit with error code if any test failed
    sys.exit(0 if failed == 0 else 1)

"""
run_tests.py
------------
Read the test registry from mhtests.yml and runs each test sequentially.

For each test, the script performs the following steps:
1. Resolves the absolute path of the test script
2. Writes this path to the minushalf.yaml file
3. Runs 'minushalf execute' within the test directory
4. Reports whether the test passes or fails

Usage:
    python run_tests.py
To add a new test, simply include a new entry in mhtests.yml. No changes to this script are required.
"""

import subprocess
import sys
import yaml  
import glob
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

import os
import shutil

def potcar_setup(registry, base_dir="test_farm/VASP"):
    """
    Copies and organizes POTCAR files for VASP test cases, based on the
    test registry loaded from the .yml file.

    For each test, reads the required atoms and:
      - Copies individual POTCAR.<atom> files to the test's minushalf_potfiles dir.
      - Concatenates them into a single POTCAR file in the test's INPUTS dir.

    Args:
        registry (list[dict]): List of test entries loaded from testconfig.yml.
        base_dir (str): Base directory where VASP test folders live.
    """
    potfiles_dir = os.path.join(base_dir, "POTFILES")

    def ensure_dir(path):
        os.makedirs(path, exist_ok=True)

    def copy_file(src, dst):
        shutil.copy2(src, dst)

    def cat_files(file_list, output_file):
        with open(output_file, 'wb') as outfile:
            for fname in file_list:
                with open(fname, 'rb') as infile:
                    outfile.write(infile.read())

    for test in registry:
        folder   = test["folder"].lstrip("/")
        atoms    = test["atoms"]
        test_dir = os.path.join(base_dir, folder.split("/")[-1])

        input_dir = os.path.join(test_dir, "INPUTS")
        mh_dir    = os.path.join(test_dir, "minushalf_potfiles")
        ensure_dir(input_dir)
        ensure_dir(mh_dir)

        potcar_paths = []
        for atom in atoms:
            src = os.path.join(potfiles_dir, f"POTCAR.{atom}")
            copy_file(src, os.path.join(mh_dir, f"POTCAR.{atom}"))
            potcar_paths.append(src)
        
        # If a *.dummy subfolder exists inside INPUTS, concatenate POTCARs there.
        # Otherwise, place the concatenated POTCAR directly in the INPUTS folder.
        dummy_dirs = glob.glob(os.path.join(input_dir, "*.dummy"))
        dummy_dirs = [d for d in dummy_dirs if os.path.isdir(d)]

        if dummy_dirs:
            for dummy_dir in dummy_dirs:
                cat_files(potcar_paths, os.path.join(dummy_dir, "POTCAR"))
        else:
            cat_files(potcar_paths, os.path.join(input_dir, "POTCAR"))

    print("POTCAR setup completed successfully.")

def load_registry(registry_file: Path):
    """Load and return the list of tests from testconfig.yml."""
    if not registry_file.exists():
        print(f"[ERROR] Registry file not found: {registry_file}")
        sys.exit(1)

    with open(registry_file, "r") as f:
        data = yaml.safe_load(f)

    # Returns a list of test entry dicts or returns an empty list is there is no data. list[dict]
    return data if data else []


def run_test(test: dict, tests_root: Path):
    """
    Run a single test entry from the registry.

    Args:
        test (dict): A test entry from mhtests.yml with keys:
                     name, description, folder, script
    """
    name        = test["name"]
    description = test["description"]
    folder      = tests_root / test["folder"]
    script      = folder / test["script"]

    print(f"\n{'='*60}")
    print(f"  TEST : {name}")
    print(f"  DESC : {description}")
    print(f"  DIR  : {folder}")
    print(f"  SCRIPT: {script}")
    print(f"{'='*60}\n")

    # Verify the test folder exists
    if not folder.exists():
        raise FileNotFoundError(f"[ERROR] Test folder not found: {folder}")
    # Verify the test script exists
    if not script.exists():
        raise FileNotFoundError(f"[ERROR] Test script not found: {script}")

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


def separator():
    print(f"\n{'='*60}")
# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    
     # Directory where test farm lives (the tests/test_farm folder)
    tests_root = Path(__file__).parent.resolve()

    tests_base_dir = tests_root / "test_farm"
    # Path to the test registry
    registry_file = tests_root / "testconfig.yml"

    separator()
    print(f"  MHTEST - Minushalf Test Runner")
    print(f"  Tests root : {tests_root}")
    print(f"  Registry   : {registry_file}")
    separator()

    tests = load_registry(registry_file)
    potcar_setup(tests)
    print(f"\n[INFO] Found {len(tests)} test(s) to run.\n")

    results = {}
    for test in tests:
        passed = run_test(test, tests_base_dir)
        results[test["name"]] = passed

    # ── Final summary ────────────────────────────────────────────────────────
    separator()
    print(f"  SUMMARY")
    separator()
    for name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"  {status}  {name}")

    total  = len(results)
    passed = sum(results.values())
    failed = total - passed
    print(f"\n  Total: {total}  |  Passed: {passed}  |  Failed: {failed}")
    separator()

    # Exit with error code if any test failed
    sys.exit(0 if failed == 0 else 1)

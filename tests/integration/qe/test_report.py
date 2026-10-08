"""Fast checks for CI reporting; no QE executable or network is required."""

import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import run_cases as report
from projwfc_nosym import disable_symmetry


class ReportTests(unittest.TestCase):
    def test_baseline_gap_and_unconverged_output(self):
        text = "convergence has been achieved\nhighest occupied, lowest unoccupied level (ev): 6.1282 6.7229\nJOB DONE."
        self.assertAlmostEqual(report.parse_baseline(text), 0.5947)
        with self.assertRaises(ValueError):
            report.parse_baseline(text.replace("convergence has been achieved", "convergence NOT achieved"))

    def test_gap_cut_formats_and_invalid_results(self):
        text = "Valence correction cuts:\n\t(Si,p):3.75 a.u\nGAP: 2.01eV"
        gap, cuts = report.parse_results(text)
        self.assertEqual(gap, 2.01)
        self.assertEqual(cuts[0]["element"], "Si")
        self.assertEqual(cuts[0]["cut_au"], 3.75)
        self.assertEqual(report.parse_results(text.replace("2.01", "0"))[0], 0)
        self.assertAlmostEqual(report.parse_results(text.replace("2.01", "2.45e-05"))[0], 2.45e-5)
        for invalid in ("", text.replace("2.01", "NaN"), text.replace("2.01", "-1"),
                        text.replace("3.75", "0"), "GAP: 2.01eV", text + "\nGAP: 3eV"):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                report.parse_results(invalid)

    def test_projection_policy_only_changes_lsym(self):
        text = "&PROJWFC\n prefix='pwscf'\n lsym = .true.\n/\n"
        self.assertEqual(disable_symmetry(text), text.replace(".true.", ".false."))
        with self.assertRaises(ValueError):
            disable_symmetry("&PROJWFC\n/\n")

    def test_fixture_files_and_reference_records(self):
        manifest = json.loads((report.ROOT / "cases.json").read_text())
        self.assertEqual([c["experimental_gap_ev"] for c in manifest["cases"]], [1.17, 7.672])
        for case in manifest["cases"]:
            self.assertTrue(case["files"])
            for name in case["files"]:
                self.assertGreater((report.ROOT / case["case"] / name).stat().st_size, 0)
                self.assertNotIn(Path(name).suffix, (".xlsx", ".xls"))

    def test_missing_results_are_not_pass_and_zero_gap_survives(self):
        manifest = json.loads((report.ROOT / "cases.json").read_text())
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            self.assertEqual([r["status"] for r in report.summarize(manifest, out)], ["NOT_RUN", "NOT_RUN"])
            case = manifest["cases"][0]
            row = report.new_result(case, manifest["dataset"])
            row.update(status="PASS", minushalf_gap_ev=0, mh_minus_exp_ev=-1.17, absolute_error_ev=1.17)
            (out / "si").mkdir()
            (out / "si/result.json").write_text(json.dumps(row))
            report.summarize(manifest, out)
            with (out / "summary.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(rows[0]["minushalf_gap_ev"], "0")
            self.assertEqual(rows[0]["mh_minus_exp_ev"], "-1.17")
            self.assertEqual(rows[1]["minushalf_gap_ev"], "")

    def test_command_failure_and_timeout(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            with self.assertRaises(RuntimeError):
                report.run_command([sys.executable, "-c", "raise SystemExit(3)"], work, work / "fail.log", 3)
            with self.assertRaises(TimeoutError):
                report.run_command([sys.executable, "-c", "import time; time.sleep(5)"], work, work / "timeout.log", 0.1)

    def test_case_comparison_and_failure_isolation(self):
        manifest = json.loads((report.ROOT / "cases.json").read_text())

        def simulate(command, work, log, timeout):
            if command[0] == "minushalf":
                (work / "minushalf_results.dat").write_text("Valence correction cuts:\n(Si,p):3.75 a.u\nGAP: 2.01eV")
            else:
                log.write_text("convergence has been achieved\nhighest occupied, lowest unoccupied level (ev): 6.1282 6.7229\nJOB DONE")

        with tempfile.TemporaryDirectory() as tmp, patch.object(report, "run_command", side_effect=simulate):
            row = report.run_case(manifest["cases"][0], manifest["dataset"], Path(tmp), 10)
            self.assertEqual(row["status"], "PASS")
            self.assertAlmostEqual(row["mh_minus_exp_ev"], 0.84)
            self.assertAlmostEqual(row["absolute_error_ev"], 0.84)
            with patch.object(report, "run_command", side_effect=RuntimeError("simulated launch failure")):
                failed = report.run_case(manifest["cases"][1], manifest["dataset"], Path(tmp), 10)
            self.assertEqual(failed["status"], "BASELINE_FAILED")
            self.assertEqual(failed["mh_minus_exp_ev"], "")
            rows = report.summarize(manifest, Path(tmp))
            self.assertEqual([r["status"] for r in rows], ["PASS", "BASELINE_FAILED"])


if __name__ == "__main__":
    unittest.main()

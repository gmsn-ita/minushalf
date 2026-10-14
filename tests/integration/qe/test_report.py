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
    def test_projection_weights_and_diagnostic_collection(self):
        from capture_diagnostics import projection_details
        weights = projection_details({"Si": [0, .01, .01, .02], "C": [0, .30, .30, .36]})
        self.assertAlmostEqual(weights["Si-p"], 4)
        self.assertAlmostEqual(weights["C-p"], 96)
        with self.assertRaises(ValueError):
            projection_details({"Si": [0, 0, 0, 0]})
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            (work / "minushalf.yaml").write_text("correction:\n  correction_code: vf\n  fractional_valence_treshold: 10\n")
            (work / "initial_vbm.json").write_text(json.dumps({"contributions_percent": weights}))
            (work / "minushalf_execute.log").write_text("Computed cut for symbol='C', orbital='p', fraction_value=100: cut=2.99")
            row = {"formula": "SiC"}
            report.collect_details({"case": "not-a-fixture"}, work, row)
            self.assertIn("Si-p: 4.000%", row["vbm_contributions"])
            self.assertEqual(row["applied_fractions"], "C-p: 100%")
            self.assertEqual(row["selection_threshold_percent"], 10)
            self.assertIn("How to read", "\n".join(report.detail_lines([row])))
            (work / "initial_vbm.json").unlink()
            report.collect_details({"case": "not-a-fixture"}, work, row)
            self.assertEqual(row["vbm_contributions"], "Unavailable")

    def test_capture_records_first_projection_and_returns_original(self):
        import os
        import types
        import capture_diagnostics as capture
        projection = {"Si": [0, 1, 0, 0]}
        class FakeBands:
            def vbm_projection(self):
                return projection
        def cli():
            self.assertIs(FakeBands().vbm_projection(), projection)
            projection["Si"][1] = 2
            FakeBands().vbm_projection()
        modules = {
            "minushalf.utils.band_structure": types.SimpleNamespace(BandStructure=FakeBands),
            "minushalf.minushalf": types.SimpleNamespace(minushalf=cli),
        }
        previous = Path.cwd()
        with tempfile.TemporaryDirectory() as tmp, patch.dict(sys.modules, modules):
            try:
                os.chdir(tmp)
                capture.main()
                saved = json.loads(Path("initial_vbm.json").read_text())
                self.assertEqual(saved["raw_projection"]["Si"][1], 1)
            finally:
                os.chdir(previous)

    def test_space_group_and_all_summary_columns_are_preserved(self):
        manifest = json.loads((report.ROOT / "cases.json").read_text())
        with tempfile.TemporaryDirectory() as tmp:
            rows = report.summarize(manifest, Path(tmp))
            self.assertEqual(rows[0]["space_group"], "Fd-3m (227)")
            text = (Path(tmp) / "summary.md").read_text()
            for column in ("Space group", "Absolute error", "CUTs", "Time (s)", "lsym"):
                self.assertIn(column, text)
            with (Path(tmp) / "summary.csv").open() as stream:
                self.assertIn("space_group", csv.DictReader(stream).fieldnames)
        self.assertEqual(report.new_result({"material_id": "unknown"}, manifest["dataset"])["space_group"], "Not recorded")

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
        self.assertEqual([c["experimental_gap_ev"] for c in manifest["cases"][:2]], [1.17, 7.672])
        for case in manifest["cases"]:
            self.assertTrue(case["files"])
            for name in case["files"]:
                self.assertGreater((report.ROOT / case["case"] / name).stat().st_size, 0)
                self.assertNotIn(Path(name).suffix, (".xlsx", ".xls"))

    def test_missing_results_are_not_pass_and_zero_gap_survives(self):
        manifest = json.loads((report.ROOT / "cases.json").read_text())
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            self.assertTrue(all(r["status"] == "NOT_RUN" for r in report.summarize(manifest, out)))
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
            if command[-1] == "execute":
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
            self.assertEqual([r["status"] for r in rows[:2]], ["PASS", "BASELINE_FAILED"])

    def test_default_selection_preserves_symmetry_and_explicit_diagnostics(self):
        manifest = json.loads((report.ROOT / "cases.json").read_text())
        selected = report.select_cases(manifest)
        self.assertTrue(all(c["projwfc_lsym"] for c in selected["cases"]))
        self.assertNotIn("mgo", [c["case"] for c in selected["cases"]])
        self.assertEqual(report.select_cases(manifest, ["mgo"])["cases"][0]["case"], "mgo")
        with self.assertRaises(ValueError):
            report.select_cases(manifest, ["unknown"])

    def test_selected_orbitals_need_not_cover_every_element(self):
        manifest = json.loads((report.ROOT / "cases.json").read_text())
        case = {**manifest["cases"][0], "required_corrected_elements": ["Si", "C"]}
        def simulate(command, work, log, timeout):
            if command[-1] == "execute":
                (work / "minushalf_results.dat").write_text("Valence correction cuts:\n(Si,p):3.75 a.u\nGAP: 2.01eV")
            else:
                log.write_text("convergence has been achieved\nhighest occupied, lowest unoccupied level (ev): 6.1282 6.7229\nJOB DONE")
        with tempfile.TemporaryDirectory() as tmp, patch.object(report, "run_command", side_effect=simulate):
            row = report.run_case(case, manifest["dataset"], Path(tmp), 10)
            self.assertEqual(row["status"], "PASS")
            self.assertEqual(row["minushalf_gap_ev"], 2.01)
            self.assertEqual(row["cuts"][0]["element"], "Si")
            self.assertAlmostEqual(row["mh_minus_exp_ev"], 2.01 - case["experimental_gap_ev"])
            self.assertEqual(row["error"], "")
            saved = json.loads((Path(tmp) / case["case"] / "result.json").read_text())
            self.assertEqual(saved["minushalf_gap_ev"], 2.01)


if __name__ == "__main__":
    unittest.main()

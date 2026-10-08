"""Run the fixed QE smoke cases and collect traceable experimental comparisons.

No API keys, spreadsheet reader or material downloads are used here.
An experimental discrepancy is reported, not used as a CI pass/fail threshold.
"""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"
FIELDS = [
    "case", "material_id", "formula", "functional", "status", "qe_gap_ev",
    "minushalf_gap_ev", "experimental_gap_ev", "mh_minus_exp_ev",
    "absolute_error_ev", "cuts", "elapsed_seconds", "projwfc_lsym",
    "corrected_upf_bytes", "scientific_validation", "dataset_doi",
    "experimental_reference_doi", "source_sheet", "source_row", "error",
]


def finite_number(value):
    number = float(value.replace("D", "e").replace("d", "e"))
    if not math.isfinite(number):
        raise ValueError("Non-finite value in calculation output")
    return number


def parse_baseline(text):
    if ("JOB DONE" not in text or "convergence has been achieved" not in text
            or "convergence NOT achieved" in text or "Error in routine" in text):
        raise ValueError("QE baseline did not finish with SCF convergence")
    levels = re.findall(
        rf"highest occupied,\s*lowest unoccupied level \(ev\):\s*({NUMBER})\s+({NUMBER})",
        text, re.IGNORECASE,
    )
    if not levels:
        raise ValueError("QE occupied/unoccupied levels not found")
    valence, conduction = map(finite_number, levels[-1])
    return max(0.0, conduction - valence)


def parse_results(text):
    gaps = re.findall(rf"(?im)^\s*GAP:\s*({NUMBER})\s*eV\s*$", text)
    if len(gaps) != 1:
        raise ValueError("Expected one finite GAP value in minushalf_results.dat")
    gap = finite_number(gaps[0])
    if gap < 0:
        raise ValueError("Negative MinusHalf gap")
    cuts, correction = [], None
    for line in text.splitlines():
        if line.strip() in ("Valence correction cuts:", "Conduction correction cuts:"):
            correction = line.split()[0].lower()
        match = re.fullmatch(rf"\s*\(([A-Z][a-z]?),\s*([spdf])\):\s*({NUMBER})\s*a\.u\s*", line)
        if match:
            element, orbital, value = match.groups()
            cut = finite_number(value)
            if correction is None or cut <= 0:
                raise ValueError("Invalid correction CUT")
            cuts.append(dict(correction=correction, element=element, orbital=orbital, cut_au=cut))
    if not cuts:
        raise ValueError("No correction CUTs found")
    return gap, cuts


def new_result(case, dataset):
    result = {key: case.get(key, "") for key in FIELDS}
    result.update(status="NOT_RUN", cuts=[], scientific_validation="NOT_PERFORMED",
                  dataset_doi=dataset["doi"], source_sheet=dataset["source_sheet"])
    return result


def run_command(command, work, log, timeout):
    """Terminate the entire Linux process group on timeout (including MPI)."""
    with log.open("w") as stream:
        process = subprocess.Popen(command, cwd=work, stdout=stream,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=max(0.01, timeout))
        except BaseException as exc:
            # Cleanup also applies to interruption; it never converts it to PASS.
            try:
                os.killpg(process.pid, signal.SIGTERM)
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                pass
            except ProcessLookupError:
                pass
            finally:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
            if isinstance(exc, subprocess.TimeoutExpired):
                raise TimeoutError(f"Time budget exceeded while running {command[0]}") from exc
            raise
    if code:
        raise RuntimeError(f"{command[0]} exited with code {code}; see {log.name}")


def prepare_case(case, work):
    """Copy only the listed fixture inputs to a fresh working directory."""
    import yaml

    fixture = ROOT / case["case"]
    for name in case["files"]:
        source = fixture / name
        shutil.copy2(source, work / name)
    config_path = work / "minushalf.yaml"
    config = yaml.safe_load(config_path.read_text())
    # Preserve scientific parameters; use two MPI processes for the CI allocation.
    config["qe"].update(
        pw_command=["mpirun", "-np", "2", "pw.x"],
        ld1_command=["ld1.x"],
        virtual_v2_command=["mpirun", "-np", "2", "virtual_v2.x"],
        projwfc_command=(["projwfc.x"] if case["projwfc_lsym"] else
                         [sys.executable, str(ROOT / "projwfc_nosym.py")]),
    )
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))


def audit_corrections(case, work, cuts):
    """Byte changes are diagnostic evidence, not proof of physical correctness."""
    audit = {}
    for element in sorted({cut["element"] for cut in cuts}):
        name = f"{element}.UPF"
        original = ROOT / case["case"] / name
        final = work / "minushalf_corrected_potentials" / name
        actual = final.read_bytes() if final.is_file() else None
        audit[name] = dict(status=("MISSING" if actual is None else
                                   "CHANGED" if actual != original.read_bytes() else "UNCHANGED"))
    (work / "correction_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    return "; ".join(f"{name}:{item['status']}" for name, item in audit.items())


def run_case(case, dataset, output, timeout):
    result = new_result(case, dataset)
    case_dir = output / case["case"]
    case_dir.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    stage = "PREPARATION"
    try:
        for name in ("baseline", "minushalf"):
            work = case_dir / name
            work.mkdir()
            prepare_case(case, work)
        stage = "BASELINE"
        baseline = case_dir / "baseline"
        log = baseline / "qe.out"
        run_command(["mpirun", "-np", "2", "pw.x", "-in", case["input_file"]],
                    baseline, log, min(180, timeout))
        result["qe_gap_ev"] = parse_baseline(log.read_text(errors="replace"))
        stage = "MINUSHALF"
        work = case_dir / "minushalf"
        remaining = timeout - (time.monotonic() - started)
        if remaining <= 0:
            raise TimeoutError("Case time budget exhausted after baseline")
        run_command(["minushalf", "execute"], work, work / "minushalf_execute.log", remaining)
        gap, cuts = parse_results((work / "minushalf_results.dat").read_text())
        result.update(minushalf_gap_ev=gap, cuts=cuts, status="PASS")
        difference = gap - case["experimental_gap_ev"]
        result.update(mh_minus_exp_ev=difference, absolute_error_ev=abs(difference),
                      corrected_upf_bytes=audit_corrections(case, work, cuts))
    except Exception as exc:
        result.update(status=f"{stage}_{'TIMEOUT' if isinstance(exc, TimeoutError) else 'FAILED'}",
                      error=str(exc))
    result["elapsed_seconds"] = round(time.monotonic() - started, 3)
    (case_dir / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"[{case['case']}] {result['status']}: MH={result['minushalf_gap_ev']} eV", flush=True)
    return result


def markdown_cell(value):
    return str(value).replace("|", "\\|").replace("\n", " ")


def format_number(value):
    return "—" if value in (None, "") else f"{value:.6g}"


def summarize(manifest, output):
    rows = []
    for case in manifest["cases"]:
        path = output / case["case"] / "result.json"
        row = new_result(case, manifest["dataset"])
        if path.is_file():
            try:
                row.update(json.loads(path.read_text()))
            except (ValueError, OSError) as exc:
                row.update(status="REPORT_FAILED", error=str(exc))
        rows.append(row)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({**row, "cuts": json.dumps(row["cuts"], separators=(",", ":"))})
    lines = ["## QE / MinusHalf reference comparison", "",
             "All gaps and errors are in eV; CUTs are in bohr (a.u.).", "",
             "| Material | XC | Status | QE gap | MH gap | Exp. gap [1] | MH − Exp. | Absolute error | CUTs | Time (s) | lsym |",
             "|---|---|---|---:|---:|---:|---:|---:|---|---:|---|"]
    for row in rows:
        cuts = "; ".join(f"{c['element']}-{c['orbital']} ({c['correction']}): {c['cut_au']:g}"
                         for c in row["cuts"]) or "—"
        cells = [row["formula"], row["functional"], row["status"],
                 *[format_number(row[k]) for k in ("qe_gap_ev", "minushalf_gap_ev", "experimental_gap_ev",
                                                   "mh_minus_exp_ev", "absolute_error_ev")],
                 cuts, format_number(row["elapsed_seconds"]), str(row["projwfc_lsym"]).lower()]
        lines.append("| " + " | ".join(map(markdown_cell, cells)) + " |")
    lines += ["", "PASS means technical completion: converged baseline, successful MinusHalf execution, finite nonnegative gap and positive CUTs.",
              "Experimental differences are descriptive; no experimental tolerance is used to pass or fail CI. Scientific validation is NOT_PERFORMED.",
              "These small fixed inputs are not a convergence study or a statistically representative benchmark. Si uses LDA/PZ; MgO uses PBE.", "",
              "[1] Borlido et al., J. Chem. Theory Comput. **15**, 5069–5079 (2019), "
              "[doi:10.1021/acs.jctc.9b00322](https://doi.org/10.1021/acs.jctc.9b00322), "
              "supporting data `ct9b00322_si_002.xlsx`, Sheet1, rows 378 (Si) and 301 (MgO). "
              "Both use the header reference: O. Madelung, *Semiconductors: Data Handbook* (2004), "
              "[doi:10.1007/978-3-642-18865-7](https://doi.org/10.1007/978-3-642-18865-7).", "",
              "### Diagnostics", ""]
    for row in rows:
        lines.append(f"- **{row['formula']}**: {markdown_cell(row['error'] or row['corrected_upf_bytes'] or row['status'])}")
    lines += ["", "UPF byte changes are recorded for inspection; a changed file alone does not validate the correction.", ""]
    (output / "summary.md").write_text("\n".join(lines))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("ci-artifacts/qe"))
    parser.add_argument("--timeout", type=float, default=None, help="Override the per-case time budgets, in seconds")
    parser.add_argument("--summarize", action="store_true", help="Collect existing reports without running calculations")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "cases.json").read_text())
    output = args.output.resolve()
    if args.summarize:
        summarize(manifest, output)
        return 0
    if args.timeout is not None and (not math.isfinite(args.timeout) or args.timeout <= 0):
        parser.error("--timeout must be finite and positive")
    if any((output / c["case"]).exists() for c in manifest["cases"]):
        parser.error("Case output already exists; choose a new --output")
    summarize(manifest, output)
    for case in manifest["cases"]:
        run_case(case, manifest["dataset"], output, args.timeout if args.timeout is not None else case["timeout_seconds"])
        rows = summarize(manifest, output)
    print((output / "summary.md").read_text())
    return 0 if all(row["status"] == "PASS" for row in rows) else 1


if __name__ == "__main__":
    sys.exit(main())

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
    "geometry_adjustment", "vbm_contributions", "selection_threshold_percent",
    "applied_fractions", "diagnostics_note",
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
        run_command([sys.executable, str(ROOT / "capture_diagnostics.py"), "execute"],
                    work, work / "minushalf_execute.log", remaining)
        gap, cuts = parse_results((work / "minushalf_results.dat").read_text())
        # Preserve numerical results even when correction coverage fails.
        result.update(minushalf_gap_ev=gap, cuts=cuts)
        difference = gap - case["experimental_gap_ev"]
        result.update(mh_minus_exp_ev=difference, absolute_error_ev=abs(difference),
                      corrected_upf_bytes=audit_corrections(case, work, cuts))
        stage = "CORRECTION_CHECK"
        required = set(case.get("required_corrected_elements", []))
        missing = required - {cut["element"] for cut in cuts}
        if missing:
            raise ValueError("Missing required corrections for: " + ", ".join(sorted(missing)))
        result["status"] = "PASS"
    except Exception as exc:
        result.update(status=f"{stage}_{'TIMEOUT' if isinstance(exc, TimeoutError) else 'FAILED'}",
                      error=str(exc))
    # Keep compact atomic/projection evidence even after a failed calculation.
    # Flattening into a visible directory allows Actions to upload these files.
    work = case_dir / "minushalf"
    evidence = case_dir / "correction-inputs"
    if work.is_dir():
        for source in work.rglob("*"):
            if source.is_file() and source.name in ("INP", "proj.in", "proj-used.in"):
                evidence.mkdir(exist_ok=True)
                relative = source.relative_to(work)
                name = "__".join(relative.parts)
                shutil.copy2(source, evidence / name.lstrip("."))
    collect_details(case, work, result)
    result["elapsed_seconds"] = round(time.monotonic() - started, 3)
    (case_dir / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"[{case['case']}] {result['status']}: MH={result['minushalf_gap_ev']} eV", flush=True)
    return result


def collect_details(case, work, result):
    """Missing evidence stays explicit; never infer contributions from final projections."""
    import yaml
    result.update(geometry_adjustment="Not documented", vbm_contributions="Unavailable",
                  selection_threshold_percent="", applied_fractions="Unavailable", diagnostics_note="")
    try:
        note = ROOT / case["case"] / "geometry-adjustment.json"
        if note.exists():
            metadata = json.loads(note.read_text())
            result["geometry_adjustment"] = metadata["description"]
            actual_input = (work / case["input_file"]).read_text()
            if not re.search(r"\bibrav\s*=\s*" + str(metadata["expected_ibrav"]) + r"\b", actual_input, re.I):
                result["geometry_adjustment"] = "Metadata does not match current ibrav; review fixture"
        config = yaml.safe_load((work / "minushalf.yaml").read_text())
        correction = config.get("correction", {})
        result["selection_threshold_percent"] = (
            correction.get("fractional_valence_treshold", "Default (not recorded)")
            if "vf" in str(correction.get("correction_code", "")) else "N/A (nonfractional mode)")
        evidence = work / "initial_vbm.json"
        if evidence.exists():
            data = json.loads(evidence.read_text())
            if "contributions_percent" in data:
                result["vbm_contributions"] = "; ".join(
                    f"{key}: {value:.3f}%" for key, value in sorted(data["contributions_percent"].items()))
            else:
                result["diagnostics_note"] = data.get("error", "Projection unavailable")
        log = work / "minushalf_execute.log"
        if log.exists():
            matches = re.findall(
                r"Computed cut for symbol='([^']+)', orbital='([^']+)', fraction_value=([0-9.]+)",
                log.read_text(errors="replace"))
            if matches:
                result["applied_fractions"] = "; ".join(
                    f"{element}-{orbital}: {fraction}%" for element, orbital, fraction in matches)
    except Exception as exc:
        result["diagnostics_note"] = str(exc)


def detail_lines(rows):
    lines = ["", "### Input and correction details", "",
             "| Material | Geometry adjustment | Initial VBM contributions | Selection threshold (%) | Applied fractions (%) |",
             "|---|---|---|---|---|"]
    for row in rows:
        cells = [row["formula"], *[("Unavailable" if row.get(key) in (None, "") else row[key]) for key in
                 ("geometry_adjustment", "vbm_contributions", "selection_threshold_percent", "applied_fractions")]]
        lines.append("| " + " | ".join(map(markdown_cell, cells)) + " |")
    for row in rows:
        if row.get("diagnostics_note"):
            lines.append("")
            lines.append(f"Diagnostic note ({row['formula']}): {markdown_cell(row['diagnostics_note'])}")
            lines.append("")
    lines += ["", "**How to read these fields**", "",
              "- **Geometry adjustment:** documented change to the fixture geometry representation. Not documented does not mean unchanged. The SiC hexagonal representation change is not a relaxation; fractional positions were retained.",
              "- **Initial VBM contributions:** normalized element/orbital weights from the first VBM projection read by MinusHalf, before correction. Full precision is saved in initial_vbm.json; displayed values are rounded. These are projected weights, not atomic charges.",
              "- **Selection threshold:** fractional-valence orbital selection cutoff. MinusHalf groups and rounds projection percentages before applying its strict greater-than threshold. Nonfractional mode is marked N/A.",
              "- **Applied fractions:** fraction_value reported for completed orbital CUT optimizations, after selection/normalization. This is a percentage of half an electron before amplitude scaling, not the original VBM weight or a count of removed electrons. Missing entries are not assumed to be zero.",
              "- **Unavailable:** the required evidence was not captured; older artifacts are not backfilled with estimates.",
              "- **CORRECTION_CHECK_FAILED:** numerical results were collected, but required element coverage failed; CI remains unsuccessful.", ""]
    return lines


def markdown_cell(value):
    return str(value).replace("|", "\\|").replace("\n", " ")


def format_number(value):
    return "—" if value in (None, "") else f"{value:.6g}"


def print_console_summary(rows):
    """Print aligned columns without Markdown markup in terminal logs."""
    headers = ("Material", "XC", "Status", "QE gap", "MH gap", "Exp. gap", "MH - Exp.")
    table = []
    for row in rows:
        values = [format_number(row[key]) for key in
                  ("qe_gap_ev", "minushalf_gap_ev", "experimental_gap_ev")]
        difference = row["mh_minus_exp_ev"]
        signed = "—" if difference in (None, "") else f"{difference:+.6g}"
        table.append([row["formula"], row["functional"], row["status"], *values, signed])
    widths = [max(len(str(record[i])) for record in [headers, *table])
              for i in range(len(headers))]

    def line(record):
        return "  ".join(str(value).ljust(widths[i]) if i < 3 else
                         str(value).rjust(widths[i]) for i, value in enumerate(record))

    print("\nQE / MinusHalf reference comparison (all gaps in eV)\n")
    print(line(headers))
    print("  ".join("-" * width for width in widths))
    for record in table:
        print(line(record))
    print("\nCalculation details:")
    for row in rows:
        cuts = "; ".join(f"{c['element']}-{c['orbital']} ({c['correction']}): {c['cut_au']:g} bohr"
                         for c in row["cuts"]) or "—"
        print(f"  {row['formula']}: CUT={cuts}; time={format_number(row['elapsed_seconds'])} s; "
              f"lsym={str(row['projwfc_lsym']).lower()}")
        if row["error"] or row["corrected_upf_bytes"]:
            print(f"    {row['error'] or row['corrected_upf_bytes']}")
    for row in rows:
        print(f"  {row['formula']}: initial VBM={row.get('vbm_contributions', 'Unavailable')}; applied={row.get('applied_fractions', 'Unavailable')}")
    print("\nPASS = technical completion; scientific validation NOT_PERFORMED.")
    print("Experimental differences are descriptive, not CI pass/fail criteria.")
    print("References and full details: summary.md and summary.csv\n")


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
    lines += detail_lines(rows)
    lines += ["", "PASS means technical completion: converged baseline, successful MinusHalf execution, finite nonnegative gap and positive CUTs.",
              "Experimental differences are descriptive; no experimental tolerance is used to pass or fail CI. Scientific validation is NOT_PERFORMED.",
              "These fixed inputs are not a convergence study or a statistically representative benchmark. See the XC column for each functional.", "",
              "[1] Borlido et al., J. Chem. Theory Comput. **15**, 5069–5079 (2019), "
              "[doi:10.1021/acs.jctc.9b00322](https://doi.org/10.1021/acs.jctc.9b00322), "
              "supporting data `ct9b00322_si_002.xlsx`, Sheet1 (individual rows in summary.csv). "
              "Default experimental reference: O. Madelung, *Semiconductors: Data Handbook* (2004), "
              "[doi:10.1007/978-3-642-18865-7](https://doi.org/10.1007/978-3-642-18865-7).", "",
              "### Diagnostics", ""]
    for row in rows:
        lines.append(f"- **{row['formula']}**: {markdown_cell(row['error'] or row['corrected_upf_bytes'] or row['status'])}")
    lines += ["", "UPF byte changes are recorded for inspection; a changed file alone does not validate the correction.", ""]
    (output / "summary.md").write_text("\n".join(lines))
    return rows


def select_cases(manifest, requested=None):
    """Default to symmetry-enabled CI fixtures; diagnostics require selection."""
    cases = manifest["cases"]
    if requested:
        unknown = set(requested) - {c["case"] for c in cases}
        if unknown:
            raise ValueError("Unknown cases: " + ", ".join(sorted(unknown)))
        selected = [c for c in cases if c["case"] in requested]
    else:
        selected = [c for c in cases if c.get("ci_enabled", True)]
        if any(not c["projwfc_lsym"] for c in selected):
            raise ValueError("Default CI cases must use lsym=true")
    if not selected:
        raise ValueError("No cases selected")
    return {**manifest, "cases": selected}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("ci-artifacts/qe"))
    parser.add_argument("--timeout", type=float, default=None, help="Override the per-case time budgets, in seconds")
    parser.add_argument("--summarize", action="store_true", help="Collect existing reports without running calculations")
    parser.add_argument("--case", action="append", help="Select a named case; repeat for multiple cases")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "cases.json").read_text())
    try:
        manifest = select_cases(manifest, args.case)
    except ValueError as exc:
        parser.error(str(exc))
    output = args.output.resolve()
    if args.summarize:
        print_console_summary(summarize(manifest, output))
        return 0
    if args.timeout is not None and (not math.isfinite(args.timeout) or args.timeout <= 0):
        parser.error("--timeout must be finite and positive")
    if any((output / c["case"]).exists() for c in manifest["cases"]):
        parser.error("Case output already exists; choose a new --output")
    summarize(manifest, output)
    for case in manifest["cases"]:
        run_case(case, manifest["dataset"], output, args.timeout if args.timeout is not None else case["timeout_seconds"])
        rows = summarize(manifest, output)
    print_console_summary(rows)
    return 0 if all(row["status"] == "PASS" for row in rows) else 1


if __name__ == "__main__":
    sys.exit(main())

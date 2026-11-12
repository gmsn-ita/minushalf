"""Download and prepare the 4H-SiC candidate locally, never during Actions.

Requires the existing preparation environment (mp-api, pymatgen, PyYAML).
The API key is read from MP_API_KEY or a hidden prompt and is never saved.
"""

import argparse
import getpass
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
import urllib.request
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
PSEUDOS = {
    element: f"https://pseudopotentials.quantum-espresso.org/upf_files/{element}.pbe-n-rrkjus_psl.1.0.0.UPF"
    for element in ("C", "Si")
}


def check_pseudo(data, element):
    upf = ET.fromstring(data)
    header = upf.find("PP_HEADER")
    if header is None or header.get("element", "").strip() != element:
        raise ValueError(f"Incorrect UPF element: expected {element}")
    if header.get("functional", "").strip() != "PBE" or header.get("is_paw") != "false":
        raise ValueError("Expected a non-PAW PBE potential")
    mesh = upf.find("PP_MESH/PP_R")
    if mesh is None or not mesh.text or float(mesh.text.split()[0]) <= 0:
        raise ValueError("Radial grid must start above zero")
    cutoffs = tuple(float(header.get(key, "0")) for key in ("wfc_cutoff", "rho_cutoff"))
    if not all(math.isfinite(value) and value > 0 for value in cutoffs):
        raise ValueError("Missing positive UPF cutoff recommendations")
    return cutoffs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kspacing", type=float, default=0.2,
                        help="Maximum reciprocal grid spacing in inverse angstrom (includes 2*pi)")
    args = parser.parse_args()
    if not math.isfinite(args.kspacing) or args.kspacing <= 0:
        parser.error("kspacing must be finite and positive")
    destination = ROOT / "sic"
    manifest_path = ROOT / "cases.json"
    manifest = json.loads(manifest_path.read_text())
    if destination.exists() or any(c["case"] == "sic" for c in manifest["cases"]):
        parser.error("SiC already exists; existing inputs were not overwritten")

    from mp_api.client import MPRester
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    import yaml

    key = os.environ.get("MP_API_KEY") or getpass.getpass("Materials Project API key (hidden): ")
    if not key:
        parser.error("A Materials Project API key is required")
    print("Downloading final structure: mp-11714 (4H-SiC)", flush=True)
    with MPRester(key) as client:
        structure = client.get_structure_by_material_id(
            "mp-11714", final=True, conventional_unit_cell=False)
    composition = structure.composition.get_el_amt_dict()
    symmetry = SpacegroupAnalyzer(structure, symprec=1e-3)
    if composition != {"Si": 4.0, "C": 4.0} or symmetry.get_space_group_number() != 186:
        raise ValueError("Expected eight-atom 4H-SiC, space group 186; review the downloaded structure")
    # Preserve the downloaded relaxed geometry. No symmetry refinement or relaxation.
    grid = [max(1, math.ceil(length / args.kspacing))
            for length in structure.lattice.reciprocal_lattice.abc]
    with tempfile.TemporaryDirectory(prefix="sic-prepare-") as tmp:
        work = Path(tmp)
        cutoff_pairs = []
        for element, url in PSEUDOS.items():
            print(f"Pseudo {element}: {url}", flush=True)
            with urllib.request.urlopen(url, timeout=120) as response:
                data = response.read()
            cutoff_pairs.append(check_pseudo(data, element))
            (work / f"{element}.UPF").write_bytes(data)
        # Round upward from UPF machine-readable values, not the rounded header prose.
        ecutwfc = math.ceil(max(pair[0] for pair in cutoff_pairs))
        ecutrho = math.ceil(max(pair[1] for pair in cutoff_pairs))
        lines = ["&CONTROL", " calculation='scf', prefix='pwscf', outdir='.',",
                 " pseudo_dir='.', verbosity='high',", "/", "&SYSTEM",
                 f" ibrav=0, nat=8, ntyp=2, ecutwfc={ecutwfc}, ecutrho={ecutrho}, nbnd=24,",
                 " nspin=1, occupations='fixed',", "/", "&ELECTRONS", "/", "",
                 "ATOMIC_SPECIES", "Si 28.0855 Si.UPF", "C 12.0107 C.UPF", "",
                 "CELL_PARAMETERS angstrom"]
        lines += [" ".join(f"{float(x):.16g}" for x in vector)
                  for vector in structure.lattice.matrix]
        lines += ["", "ATOMIC_POSITIONS crystal"]
        lines += [site.specie.symbol + " " + " ".join(f"{float(x):.16g}" for x in site.frac_coords)
                  for site in structure]
        lines += ["", "K_POINTS automatic", " ".join(map(str, grid)) + " 0 0 0", ""]
        (work / "qe.in").write_text("\n".join(lines))
        config = {
            "software": "QE",
            "qe": {"pw_command": ["mpirun", "-np", "2", "pw.x"],
                   "projwfc_command": ["projwfc.x"], "ld1_command": ["ld1.x"],
                   "virtual_v2_command": ["mpirun", "-np", "2", "virtual_v2.x"],
                   "input_file": "./qe.in"},
            "atomic_program": {"exchange_correlation_code": "pb", "calculation_code": "ae", "max_iterations": 100},
            "correction": {"correction_code": "vf", "potfiles_folder": "./",
                           "amplitude": 1.0, "indirect": ".true.", "fractional_valence_treshold": 10},
        }
        (work / "minushalf.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
        (work / ".gitattributes").write_text("*.UPF -whitespace\n")
        case = {
            "case": "sic", "material_id": "mp-11714", "formula": "4H-SiC",
            "input_file": "qe.in", "functional": "PBE", "projwfc_lsym": True,
            "ci_enabled": True, "experimental_gap_ev": 3.2,
            "experimental_reference_doi": "10.1007/978-3-642-18865-7",
            "reference_source": "Madelung handbook, as specified in the workbook header; no individual DOI in this row.",
            "source_row": 383, "reference_icsd_id": "24170",
            "scientific_validation": "NOT_PERFORMED", "required_corrected_elements": ["Si", "C"],
            "files": ["Si.UPF", "C.UPF", "qe.in", "minushalf.yaml"], "timeout_seconds": 1800,
        }
        manifest["cases"].append(case)
        shutil.copytree(work, destination)
        try:
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        except BaseException:
            shutil.rmtree(destination)
            raise
    print(f"Prepared candidate: {destination.name}; k grid {grid}; cutoffs {ecutwfc}/{ecutrho} Ry")
    print("vf uses computed VBM contributions and a 10% orbital threshold; fractions are NOT hardcoded.")
    print("Run the candidate and inspect selected orbitals, fractions and atomic inputs before accepting it.")
    print("No calculations executed; lsym=true has not yet been validated for this fixture.")


if __name__ == "__main__":
    main()

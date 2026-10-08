# Small QE integration cases

These fixed Si and MgO cases exercise the QE/MinusHalf workflow after a code
change. Keeping the inputs in the repository makes each run repeatable and
avoids Materials Project API keys or pseudopotential downloads during CI.
The container build and Python installation still download their dependencies.

## Experimental references

Only the following two entries were transcribed from the supporting workbook;
the XLSX is neither included in the repository nor required by the tests.

| Case | Reference MP ID | Reference ICSD ID | Experimental gap (eV) | Sheet1 row | Source |
|---|---|---|---:|---:|---|
| Si | mp-149 | 51688 | 1.17 | 378 | [1], with experimental reference [2] in the workbook header |
| MgO | mp-1265 | 9863 | 7.672 | 301 | [1], with experimental reference [2] in the workbook header |

1. P. Borlido, T. Aull, A. W. Huran, F. Tran, M. A. L. Marques and S. Botti,
   *Large-Scale Benchmark of Exchange–Correlation Functionals for the
   Determination of Electronic Band Gaps of Solids*, Journal of Chemical
   Theory and Computation **15**, 5069–5079 (2019).
   DOI: [10.1021/acs.jctc.9b00322](https://doi.org/10.1021/acs.jctc.9b00322).
   Supporting workbook: `ct9b00322_si_002.xlsx`, `Sheet1`.
2. O. Madelung, *Semiconductors: Data Handbook* (2004).
   DOI: [10.1007/978-3-642-18865-7](https://doi.org/10.1007/978-3-642-18865-7).
   This is the workbook's default reference; neither selected row supplies
   a separate experimental DOI. These values are transcribed from [1], not
   independently re-extracted from the handbook.

`cases.json` contains the two reference values, citations, case settings and
input filenames. Download metadata and checksums are kept outside Git.

## Inputs and pseudopotentials

### Silicon

The existing diamond-Si input and `Si.UPF` are retained byte for byte from
commit `e3f758e`. This is the validated CI fixture, not a new download of the
Materials Project entry or an exact reproduction of the benchmark structure.
Its MP ID above identifies the corresponding experimental reference entry.

- Primitive two-atom cell, `ibrav=2`, `celldm(1)=10.2076` bohr.
- `8 8 8 0 0 0` k-point grid; `ecutwfc=50` Ry; eight bands.
- The committed `Si.UPF` header specifies a norm-conserving LDA/PZ potential.
- Atomic correction uses `pz`, valence correction, amplitude 1 and indirect gap.
- `projwfc` uses `lsym=.true.`.

The inherited Si UPF does not name its original author/download URL. Its
provenance is the existing repository fixture; no PSlibrary
or SSSP source is claimed for it.

### Magnesium oxide

The rocksalt MgO input and pseudopotentials are copied unchanged from the
locally downloaded and prepared Materials Project case `mp-1265`. The cell and
fractional positions were checked against the downloaded structure before
including the input. The k-point grid remains `13 13 13 0 0 0`.

- Two atoms, explicit cell in angstrom; twelve bands and fixed occupations.
- `ecutwfc=59` Ry; `ecutrho=406.723404255` Ry, as in the pilot.
- PBE ultrasoft pseudopotentials from the official QE PSlibrary download site:
  - `Mg.UPF`: [Mg.pbe-spnl-rrkjus_psl.1.0.0.UPF](https://pseudopotentials.quantum-espresso.org/upf_files/Mg.pbe-spnl-rrkjus_psl.1.0.0.UPF).
  - `O.UPF`: [O.pbe-n-rrkjus_psl.1.0.0.UPF](https://pseudopotentials.quantum-espresso.org/upf_files/O.pbe-n-rrkjus_psl.1.0.0.UPF).
- These two files were downloaded once and committed unchanged under shorter
  names. Their original headers, credits and generation data are preserved.
- Atomic correction uses `pb` (PBE), valence correction, amplitude 1 and indirect gap.
- `projwfc_nosym.py` explicitly uses `lsym=.false.` for MgO, following the
  previously tested workaround for the `d_matrix` non-orthogonality error.
  It saves the actual input as `proj-used.in` in each calculation directory.
  This affects the projection step; it does not set `nosym` in the SCF input.

Neither k-point mesh has been established as converged here. The symmetry
workaround and the difference between the two XC choices must be considered
before using these numbers as scientific results.

The runner copies fixtures into separate baseline and MinusHalf directories,
leaving the repository inputs unchanged. In the working YAML it sets the
executable commands for the container: two MPI processes for `pw.x` and
`virtual_v2.x`, and the explicit per-case projection policy. Scientific
parameters from the fixture YAML are retained. The working YAML is saved
with the results.

## Results and checks

The Actions run summary displays both materials. Download its artifact for
`ci-artifacts/qe/summary.csv`, per-case `result.json`, logs and input/version
records. A case that fails does not prevent the next case from being attempted;
the overall integration step still fails when any case fails. Unit tests must
pass before the material runs begin.

The columns include:

- `qe_gap_ev`: unoccupied minus occupied level from the converged baseline QE output.
- `minushalf_gap_ev`: the gap reported in `minushalf_results.dat`.
- `experimental_gap_ev`: the transcribed reference value.
- `mh_minus_exp_ev`: MinusHalf minus experiment; negative means underestimation.
- `absolute_error_ev`: the absolute value of that difference.
- `cuts`: element, orbital, valence/conduction correction and radius in bohr.
- `elapsed_seconds`: total case time, including the baseline and MinusHalf run.
- `projwfc_lsym`: the explicitly configured projection policy.
- `corrected_upf_bytes`: whether the exported corrected UPFs changed compared
  with the inputs, or were missing. The status is saved in `correction_audit.json`.
- Dataset DOI, experimental reference DOI, source sheet/row and any error.

Missing results are blank, not zero. A real zero gap is retained. No experimental
tolerance controls PASS/FAIL. PASS requires a converged, completed baseline,
a zero exit status from MinusHalf, one finite nonnegative reported gap and
at least one finite positive CUT. It does not prove that the correction is
physically correct. `scientific_validation` remains `NOT_PERFORMED`.

Before adding physical acceptance thresholds, check the structure/polymorph,
experimental conditions and type of gap, convergence and the actual corrected
potential. The grid-sampled baseline gap is not a full band-extrema search.
Do not tune these inputs merely to make the experimental error small.

## Running and extending

After installing MinusHalf and QE in the existing test container:

```bash
python tests/integration/qe/run_cases.py
python tests/integration/qe/run_cases.py --summarize
```

For another local run, supply a fresh `--output` directory. Existing case
directories are not overwritten. Si has a 900-second budget and MgO has
1800 seconds, including their baseline (capped at 180 seconds). `--timeout`
overrides the total per-case budgets. The Actions job has a 60-minute limit
including the build and installation. These are CI budgets, not cluster settings.

To add a material, add a fixture directory and a `cases.json` entry with its
reference, input filenames and time budget, document its origin and confirm its runtime.
The collector uses the manifest to include all cases. SiC is a candidate for
a later addition after confirming that its polytype matches the experimental
reference and testing the fractional correction configuration with the group.

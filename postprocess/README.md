# Sputtering Yield & Molecule Analysis

Post-processing pipeline for LAMMPS molecular dynamics simulations of ion bombardment / sputtering. Given a set of per-impact simulation folders (each containing LAMMPS dump files and ingress/egress count data), it computes:

- **Sputtering yields** (physical, chemical, total) with statistical error estimates, per ion energy and per simulation directory.
- **Sputtered molecule identification**: detects diatomic and polyatomic target-target, target-ion, and ion-ion species leaving the surface, based on interatomic distance and potential-energy thresholds.
- **Rovibrational analysis**: for identified diatomic species, decomposes kinetic energy into center-of-mass translational, rotational, and vibrational components, and estimates rotational (J) and vibrational (n) quantum numbers using rigid-rotor and harmonic-oscillator approximations, calibrated against Tersoff-potential-derived bond parameters (`vibrational_data` in the script).

## Requirements

```
pip install -r requirements.txt
```

Needs Python 3.10+ (uses PEP 572 walrus operator) and pandas >= 2.2 (uses `include_groups` in `groupby().apply()`).

## Usage

```
python3 analysis.py
```

The script prompts interactively for:
1. Target element symbol (e.g. `W`) and ion element symbol (e.g. `D`).
2. One or more `(label, directory path)` pairs — each directory should contain energy-named subfolders (e.g. `D1000`, `D3000`, ...), each holding numbered simulation-run subfolders (e.g. `1`, `2`, ...) with LAMMPS dump `.txt` files and per-run `.csv` ingress/egress data.

Output: `sputtering_yields_<ion><target>.json` (yields per directory/energy), plus per-energy-folder CSVs (`Target_molecules.csv`, `diatomic_target.csv`, `polyatomic_target.csv`, `ion_molecules.csv`, `ion_single.csv`, `final_rovib_target.csv`, `final_rovib_ion.csv`, `event.csv`, `sputtered.data`) and a `<ion><target>_<label>_sputtered_species.csv` summary per directory.

## Notes

- Uses `ProcessPoolExecutor`/`ThreadPoolExecutor` at several stages for parallelism across simulation folders and files; scales with `os.cpu_count()`.
- `vibrational_data` currently covers Be/W targets against H/D/T ions — extend this table for other target/ion combinations.

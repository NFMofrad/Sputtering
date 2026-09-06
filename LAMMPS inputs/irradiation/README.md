# irradiation

Fires a single D ion at an equilibrated W surface (from `../surface_gen/`) and
runs the impact until the ion/sputtered atoms exit a control region or a time
limit is hit. Meant to be run once per row of the `.tsv` produced by
`generate_input.py` — i.e. thousands of independent single-impact runs per
energy/angle set.

## Files

| File | Purpose |
|---|---|
| `generate_input.py` | Generates initial ion position/velocity conditions for a batch of impacts at a given energy and angle, written to a `.tsv`. |
| `irrad.lmp` | LAMMPS script for one ion-impact simulation. |
| `elstop-DW.txt` | Electronic stopping power table for D in W (used by `fix electron/stopping`). |

## Requirements

- LAMMPS with `MANYBODY` (`pair_style tersoff/zbl`), `MISC` (`fix
  electron/stopping`, `fix halt`, `fix dt/reset`), and `EXTRA-FIX` as needed
  for your LAMMPS version.
- Python 3.10+ with:
  ```
  numpy
  pandas
  periodictable
  ```
  (install via the top-level `LAMMPS inputs/requirements.txt`, or `pip
  install numpy pandas periodictable`).
- `WD.tersoff.zbl` from `../../potentials/`, copied into the run directory.
- `W1D1_open.eq_300.restart`, produced by `../surface_gen/surface_gen.lmp`
  (see that folder's README) — must be present in the run directory.

## Step 1 — generate impact conditions

```bash
python3 generate_input.py --energy 2000 --theta 0 --ion D --n-phi-samples 12000
```

Key options:

| Flag | Meaning | Default |
|---|---|---|
| `--energy` | Ion energy, eV | `2000` |
| `--theta` | Incidence angle, deg | `0` |
| `--ion` | Ion element symbol (mass looked up via `periodictable`) | `D` |
| `--dist` | Initial height above surface, Å | `8` |
| `--n-phi-samples` | Number of independent impact runs to generate | `12000` |
| `--seed` | RNG seed | `190820222` |
| `--output` | Output path | `Input_<energy>eV_<theta>deg.tsv` |

Output: a tab-separated file with one row per run — `run_id, theta, phi, x0,
y0, z, v, vx, vy, vz, shift_x, shift_y`.

## Step 2 — run LAMMPS per row (not included — you need to write this)

**There is currently no driver script in the repo that turns `.tsv` rows into
LAMMPS runs.** `irrad.lmp` expects these variables to be passed on the
command line, one set per impact:

| Variable | Source column in the `.tsv` |
|---|---|
| `T` | fixed per batch, e.g. `300` (must match the restart file's temperature) |
| `vx`, `vy`, `vz` | `vx`, `vy`, `vz` |
| `shiftx`, `shifty` | `shift_x`, `shift_y` |
| `i` | `run_id` (used to name `control.${i}.csv`) |
| `rid` | `run_id` (logged in the thermo CSV) |

A minimal loop (e.g. `run_batch.sh`) would look like:

```bash
#!/usr/bin/env bash
i=4
TSV_FILE="$PWD/Input_2000eV_0deg.tsv"

# Read position and velocity from TSV based on array task ID
line=$(sed -n "${i}p" "$TSV_FILE")
IFS=$'\t' read -r rid theta phi x y z v vx vy vz shiftx shifty <<< "$line"

# Run simulation
mkdir $rid
cd $rid
lmp -in ../in.lmp -log Irrad.$rid.log -var i $rid -var rid $rid -var T 300 -var x $x -var y $y -var z $z -var vx $vx -var vy $vy -var vz $vz -var shiftx $shiftx -var shifty $shifty
cd ..
```

On an HPC cluster you'd more likely submit this as a job array (one task per
`run_id`) — adapt to your scheduler.

## Output (per run)

- `control.${i}.txt` — per-run LAMMPS dump of atoms in the sputtering control
  region (id, element, type, position, velocity, KE, PE).
- `control.lammpstrj` — full trajectory dump (commented for production runs —
  **remove/comment the `DUMP_GROUPS` dump line in `irrad.lmp` before large
  batch runs**, or it will generate large amounts of data per impact).
- `control.${i}.csv` — thermo/count time series (temp, pressure, energies,
  atom counts in control regions, electronic stopping loss).

These per-run outputs (organized as `<label>/<energy>/<run_id>/`) are exactly
what `../../postprocess/analysis.py` expects as input.

## Notes

- `irrad.lmp` is hardcoded to W/D and `T=300` via the restart filename —
  update `read_restart` and `pair_coeff` for other targets/temperatures.
- Exit conditions: halts when no target atoms remain in the control regions
  and simulated time > 3 ps, or unconditionally after 10 ps (`FIX_HALT_TIME`).
- Uses an adaptive timestep (`fix dt/reset`).

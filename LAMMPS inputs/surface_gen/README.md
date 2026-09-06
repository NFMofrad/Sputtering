# surface_gen

Builds and equilibrates a W(110) slab with a free surface, then deposits D atoms
onto it and relaxes the result. This is the **first stage** of the pipeline —
its output restart file is the required input for `LAMMPS inputs/irradiation/`.

## Files

| File | Purpose |
|---|---|
| `surface_gen.lmp` | LAMMPS script: equilibrate bulk W, open a surface, deposit D, relax. |
| `W.data` | Initial W(110) atomic configuration (31,360 atoms), generated with [Atomsk](https://atomsk.univ-lille.fr/). |

## Requirements

- LAMMPS (built with the `MANYBODY` package for `pair_style tersoff/zbl`, and
  the `MISC` package for `fix deposit` / `fix halt` / `fix dt/reset`).
- The potential file `WD.tersoff.zbl` from `../../potentials/`, copied or
  symlinked into this directory before running.

## Variables you must define

`surface_gen.lmp` does not set any of these itself — pass them on the command
line with `-var`:

| Variable | Meaning | Example |
|---|---|---|
| `T` | Target/bath temperature, K | `300` |
| `P` | Target pressure for the barostat, bar | `0` |
| `Time` | Equilibration duration, ps (converted to steps via `ceil(Time/dt)`) | `10` |
| `dt` | Timestep, ps | `0.001` |
| `seed` | RNG seed for initial velocities and D deposition | `12345` |

## Running

```bash
lmp -in surface_gen.lmp -log eq.log -var T 300 -var P 0.0 -var dt 0.001 -var Time 30 -var seed $RANDOM
```

## Output

- `W1D1_open.eq_${T}.restart` — **this is the file `irradiation/irrad.lmp`
  reads at the top of its script** (currently hardcoded to
  `W1D1_open.eq_300.restart`, i.e. `T=300`). If you run a different `T`,
  either rename the restart file or edit the `read_restart` line in
  `irrad.lmp` to match.
- `W1D1_open.eq_${T}.data` — plain-text equivalent of the restart file.
- `W.minimized.lammpstrj`, `W_groups.eq_${T}.lammpstrj`,
  `W_final.eq_${T}.lammpstrj` — trajectory dumps at each stage (minimization,
  post-deposition groups, final relaxation).
- `W.minimized.csv` — thermo log (time, T, P, energies, box dims) during
  minimization/equilibration.

## Notes

- `thermo_modify lost warn` during deposition means atoms leaving the box are
  tolerated (not fatal) — expected behavior during the deposition run.

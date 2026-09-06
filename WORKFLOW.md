# End-to-end workflow

This repo runs as a three-stage pipeline. Each stage's output is the next
stage's required input — none of this is automatic, so files have to be
copied between directories by hand (or by a wrapper script you write).

```
LAMMPS inputs/surface_gen/          LAMMPS inputs/irradiation/            postprocess/
┌───────────────────────┐          ┌───────────────────────┐          ┌───────────────────────┐
│ W.data                │          │ generate_input.py      │          │ analysis.py           │
│ surface_gen.lmp        │          │ irrad.lmp               │          │                       │
│  + WD.tersoff.zbl  ───┼───►      │  + WD.tersoff.zbl  ◄───┼── needs  │  reads dump + csv     │
│                        │ restart  │  + W1D1_open...restart │  restart │  from irradiation/    │
│                        │  file    │  + elstop-DW.txt       │  file    │  run folders           │
└───────────────────────┘          └───────────────────────┘          └───────────────────────┘
   equilibrate slab,                  one process per ion impact,          yields, molecule ID,
   deposit D, relax                   thousands of runs per (E, θ)         rovibrational analysis
```

## 1. Build and equilibrate the surface

```bash
cd "LAMMPS inputs/surface_gen"
cp ../../potentials/WD.tersoff.zbl .
lmp -in surface_gen.lmp -var T 300 -var P 0 -var Time 10 -var dt 0.001 -var seed 12345
```
Produces `W1D1_open.eq_300.restart` (see `surface_gen/README.md`).

## 2. Generate impact conditions and run irradiation

```bash
cd "../irradiation"
python3 generate_input.py --energy 2000 --theta 0 --ion D --n-phi-samples 12000
cp ../../potentials/WD.tersoff.zbl ../surface_gen/W1D1_open.eq_300.restart .
# then run one LAMMPS invocation per row of the generated .tsv (see irradiation/README.md)
```
Organize outputs as `<label>/<energy>/<run_id>/` so `postprocess/analysis.py`
can find them.

## 3. Postprocess

```bash
cd ../../postprocess
pip install -r requirements.txt
python3 analysis.py
# interactively enter: target element, ion element, and (label, directory) pairs
```
Produces sputtering yields, molecule/rovibrational CSVs, and a JSON summary.

## What's not (yet) in the repo

- A driver/submission script connecting step 2's `.tsv` output to individual
  `irrad.lmp` runs (example given in `irradiation/README.md` — worth turning
  into a real script or Slurm/PBS array job).
- Be/D input files (only W/D `surface_gen`/`irradiation` inputs exist, though
  `potentials/BeD.tersoff.zbl` is present).
- Automated copying of potential/restart files between stages.

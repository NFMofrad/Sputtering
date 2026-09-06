# potentials

Tersoff/ZBL interatomic potential files for use with LAMMPS `pair_style
tersoff/zbl`, consumed by both `LAMMPS inputs/surface_gen/` and
`LAMMPS inputs/irradiation/`.

| File | System | Source / citation (from file header) |
|---|---|---|
| `WD.tersoff.zbl` | W–D | N. Juslin, P. Erhart, P. Träskelin, J. Nord, K. O. E. Henriksson, K. Nordlund, E. Salonen, K. Albe — DOI: 10.1063/1.2149492 |
| `BeD.tersoff.zbl` | Be–D | C Björkas et al 2009 J. Phys.: Condens. Matter 21 445002 DOI: 10.1088/0953-8984/21/44/445002 |

## Usage

```
pair_style      tersoff/zbl
pair_coeff      * * WD.tersoff.zbl W D      # for W/D systems
pair_coeff      * * BeD.tersoff.zbl Be D    # for Be/D systems
```

## Notes

- `BeD.tersoff.zbl` is present but **not yet wired into any input script** in
  this repo — `surface_gen.lmp` and `irrad.lmp` are both currently W/D only.
  To run Be systems you'd need a Be `.data` starting structure and a
  corresponding `elstop-DBe.txt` electronic-stopping table.

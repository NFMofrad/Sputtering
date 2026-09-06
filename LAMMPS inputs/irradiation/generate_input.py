#!/usr/bin/env python3
"""Generate initial position/velocity conditions for LAMMPS ion-bombardment runs."""

import argparse

import numpy as np
import pandas as pd
import periodictable


def compute_kinetic(energy, mass, dist, theta, phi, impact_x=0.0, impact_y=0.0, impact_z=0.0):
    z = dist + impact_z
    r = dist / np.cos(np.radians(theta))
    x = impact_x + (r * np.sin(np.radians(theta)) * np.cos(np.radians(phi)))
    y = impact_y + (r * np.sin(np.radians(theta)) * np.sin(np.radians(phi)))

    v = np.sqrt(2 * energy / mass / 931494102.42) * 299792458.0 * 1e-2
    vx = v * np.sin(np.pi - np.radians(theta)) * np.cos(np.pi + np.radians(phi))
    vy = v * np.sin(np.pi - np.radians(theta)) * np.sin(np.pi + np.radians(phi))
    vz = v * np.cos(np.pi - np.radians(theta))
    return x, y, z, v, vx, vy, vz


def parse_args():
    parser = argparse.ArgumentParser(description="Generate initial ion-bombardment conditions and write a .tsv input file.")
    parser.add_argument("--energy", type=float, default=2000, help="ion energy, eV (default: %(default)s)")
    parser.add_argument("--theta", type=float, default=0, help="ion incidence angle, deg (default: %(default)s)")
    parser.add_argument("--ion", type=str, default="D", help="ion element symbol, e.g. 'D' (atomic mass is looked up automatically) (default: %(default)s)")
    parser.add_argument("--dist", type=float, default=8, help="initial distance above the surface, Å (default: %(default)s)")
    parser.add_argument("--impact-x", type=float, default=0.0, help="impact point x offset, Å (default: %(default)s)")
    parser.add_argument("--impact-y", type=float, default=0.0, help="impact point y offset, Å (default: %(default)s)")
    parser.add_argument("--impact-z", type=float, default=0.0, help="impact point z offset, Å (default: %(default)s)")
    parser.add_argument("--n-phi-samples", type=int, default=12000, help="number of simulation runs / uniform random phi samples in [0, 180) deg (default: %(default)s)")
    parser.add_argument("--seed", type=int, default=190820222, help="random seed for phi sampling and position shifts (default: %(default)s)")
    parser.add_argument("--output", type=str, default=None, help="output .tsv path (default: Input_<energy>eV_<theta>deg.tsv)")
    return parser.parse_args()


def main():
    args = parse_args()

    ion = args.ion.strip().capitalize()
    try:
        mass = getattr(periodictable, ion).mass
    except AttributeError:
        raise SystemExit(f"Error: '{args.ion}' is not a valid element symbol.")
    print(f"Mass of {ion}: {mass} u")

    # Sample ion impact angle phi uniformly; theta is fixed for this input set
    ion_data = {
        "theta": args.theta,
        "phi": np.empty(args.n_phi_samples),
    }

    # build dataframe with combinations of theta and number of phi samples
    db = pd.DataFrame(
        np.array(np.meshgrid(*ion_data.values())).T.reshape(-1, len(ion_data.keys())),
        columns=ion_data.keys(),
        dtype=float,
    )

    # fill in values of phi from uniform distribution
    rng = np.random.default_rng(seed=args.seed)
    db["phi"] = rng.uniform(0, 180, size=len(db))

    # set unique simulation IDs
    db["run_id"] = range(1, len(db) + 1)

    # fill in values of position and velocity, given theta and phi
    db[["x0", "y0", "z", "v", "vx", "vy", "vz"]] = db.apply(
        lambda row: compute_kinetic(
            energy=args.energy,
            mass=mass,
            dist=args.dist,
            theta=row.theta,
            phi=row.phi,
            impact_x=args.impact_x,
            impact_y=args.impact_y,
            impact_z=args.impact_z,
        ),
        axis=1,
        result_type="expand",
    )

    db["shift_x"] = rng.uniform(-0.5, 0.5, size=len(db))
    db["shift_y"] = rng.uniform(-0.5, 0.5, size=len(db))

    # clean up and output dataframe
    db.sort_values(by=["theta", "run_id"], inplace=True)
    db.set_index("run_id", inplace=True)

    comments = {
        "run id": "-",
        "ion theta": "deg",
        "ion phi": "deg",
        "x position": "Å",
        "y position": "Å",
        "z position": "Å",
        "velocity": "Å/ps",
        "x velocity": "Å/ps",
        "y velocity": "Å/ps",
        "z velocity": "Å/ps",
        "random shift x": "-",
        "random shift y": "-",
    }

    output_path = args.output or f"Input_{args.energy:g}eV_{args.theta:g}deg.tsv"
    with open(output_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("# " + "\t".join(comments.keys()) + "\n")
        f.write("# " + "\t".join(comments.values()) + "\n")
        db.to_csv(f, sep="\t", float_format="%.9f")

    with pd.option_context("display.max_rows", 50):
        print(db)

    print(f"There are {len(db.index)} simulations in this set (written to {output_path})")


if __name__ == "__main__":
    main()

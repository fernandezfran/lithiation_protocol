#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lithiation of amorphous silicon using DFTB+.

We start with an amorphized silicon structure and follow the next protocol:

1. add a Li atom at the center of the largest spherical void,
2. increase the volume and scale the coordinates,
3. perform a local LBFGS minimization,
4. simulate 10ps with Berendsen NPT,
5. select the frame with minimum absolute pressure,
6. with x defined as number of Li atoms per Si atoms if x < 3.75 goto point 1
else finish.

To make the process faster, 4 atoms of lithium are added at a time and expanding
in each one of them, this was checked with DFT which does not alter the results.
"""
import argparse

import analysis

import lithiation


def main():
    parser = argparse.ArgumentParser(
        description="Lithiate an amorphous structure, "
        "by default from the beginning but can also be restarted from a given "
        "structure. You have also different options to analyze the lithiation "
        " once it was performed."
    )

    parser.add_argument(
        "--restart-from",
        help="restart from a given structure RESTART_FROM, e.g. Li55Si64",
    )

    parser.add_argument(
        "--fvc",
        action="store_true",
        help="fractional volume change calculation",
    )

    parser.add_argument(
        "--rdf",
        nargs="?",
        const=5,
        type=int,
        help="rdf specifying how many structures to skip for the plot, by default 5",
    )

    parser.add_argument(
        "--central",
        nargs="?",
        const="Si",
        type=str,
        help="central atom type for rdf calculation, by default Si",
    )

    parser.add_argument(
        "--interact",
        nargs="?",
        const="Si",
        type=str,
        help="interact atom type for rdf calculation, by default Si",
    )

    parser.add_argument(
        "-s",
        "--save",
        action="store_true",
        help="to save the png figures if created",
    )

    args = parser.parse_args()

    full = True
    for v in args.__dict__.values():
        if v:
            full = False

    if full or args.restart_from:
        lithiation.Lithiation(structure=args.restart_from, nsteps=4).run()

    if args.fvc:
        analysis.fvc_plot(save=args.save)

    if args.rdf:
        if args.central is None:
            args.central = "Si"
        if args.interact is None:
            args.interact = "Si"
        analysis.rdf_plot(
            each=args.rdf, central=args.central, interact=args.interact, save=args.save
        )


if __name__ == "__main__":
    main()

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
import os
import subprocess

import exma
import numpy as np

import lithiation


def main():
    parser = argparse.ArgumentParser(
        description="Lithiate an amorphous structure, "
        "by default from the beginning but can also be restarted from a given "
        "structure"
    )
    parser.add_argument(
        "--restart-from",
        help="restart from a given structure RESTART_FROM, e.g. Li55Si64",
    )
    args = parser.parse_args()

    lithiation.Lithiation(structure=args.restart_from, nsteps=4).run()


if __name__ == "__main__":
    main()

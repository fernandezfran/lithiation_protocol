#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Full lithiation of amorphous silicon using DFTB+.

We start with an amorphized silicon structure and follow the next protocol:

1. add a Li atom at the center of the largest spherical void,
2. increase the volume and scale the coordinates,
3. perform a local LBFGS minimization,
4. simulate 10ps with Berendsen NPT,
5. select the frame with minimum absolute pressure,
6. with x defined as number of Li atoms per Si atoms if x < 3.75 goto point 1
else finish.

To make the process faster, 3 atoms of lithium are added at a time and expanding
in each one of them, this was checked with DFT which does not alter the results.
"""
import os
import subprocess

import exma
import numpy as np

from lithiation_step import lithiation_step
from io_dftb_plus import read_md_out, write_gen_format


def main():
    # get the last frame of a trajectory and wrap the frame inside the box
    frames = exma.read_xyz("a-Si64.xyz")
    frame = frames[-1]
    frame.box = np.full(3, 10.937456)
    frame._wrap()

    # initial number of Li and Si atoms
    nli = 0
    nsi = frame.natoms

    # add a Li atom and expand frame
    dmax, frame = lithiation_step(frame, 10.937456)
    nli += 1

    x = nli / nsi
    while x < 3.75:
        # write genFormat file for dftb+
        write_gen_format(frame, "LixSi64.gen")

        # run LBFGS minimization and Berendsen NPT
        subprocess.run(["bash", "run.sh"])

        # get frame with minimum pressure
        df = read_md_out()
        min_press_frame = np.argmin(np.abs(df["press"].values))

        frames = exma.read_xyz("LixSi64.xyz")
        frame = frames[min_press_frame]
        frame.box = np.array(
            [df[kbox][min_press_frame] for kbox in ("xbox", "ybox", "zbox")]
        )
        frame._wrap()

        # save files of interest
        os.system(f"mv md.out dftb+npt/md.Li{nli}Si64.out")
        os.system(f"mv LixSi64.xyz dftb+npt/Li{nli}Si64.xyz")

        # add 3 Li atoms and expand frame
        for i in range(3):
            dmax, frame = lithiation_step(frame, frame.box[0])
            nli += 1
        x = nli / nsi


if __name__ == "__main__":
    main()

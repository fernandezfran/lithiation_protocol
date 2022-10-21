#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Full lithiation of amorphous silicon using DFTB+.

We start with an amorphous Si structure and follow the protocol:
1. add a Li atom at the center of the largest spherical void,
2. increase the volume and scale the coordinates,
3. simulate 5ps with Berendsen NPT,
4. select the frame with minimum absolute pressure,
5. goto point 1 if x < 3.75 else finish.

To make the process faster, 3 atoms of lithium are added at a time and expanding
in each one of them, this was checked with DFT which does not alter the results.
"""
import exma
import numpy as np

from lithiation_step import lithiation_step


def main():
    box = 10.937456

    # get the last frame of a trajectory
    frames = exma.read_xyz("a-Si64.xyz")
    frame = frames[-1]

    # wrap the frame inside the box
    frame.box = np.full(3, box)
    frame._wrap()

    dmax, frame = lithiation_step(frame, box, 0.5)


if __name__ == "__main__":
    main()

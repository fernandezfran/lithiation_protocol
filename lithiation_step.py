#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Slightly similar Chevrier and Dahn protocol for lithiation of a-Si.

I only consider the first two steps:
1. add a Li atom at the center of the largest spherical void,
2. increase the volume and scale the coordinates,
the rest (volume minimization, energy, etc.) are replaced by a NPT simulation.

https://doi.org/10.1149/1.3111037 is the paper where the protocol is established.
"""
import exma
import numpy as np


def lithiation_step(frame, box, dx):
    """A single step of the lithiation protocol.

    Add a Li atom at the center of the largest spherical void in a structure
    and expand the volume.

    To find the largest spherical void, a 3d grid of points separated by `dx` in
    each direction is generated and the distance from each of them to the
    nearest atom is found, then the largest one is selected.

    Parameters
    ----------
    frame : exma.core.AtomicSystem
        frame of the structure to lithiate and expand

    box : float
        the box size

    dx : float
        separation between 3d grid points in each direction

    Returns
    -------
    tuple
        with a float corresponding with the radius of the sphere and a frame
        (exma.core.AtomicSystem) with the new structure
    """
    # generate the grid of points to center spherical voids
    grid1d = np.linspace(0, box, num=np.intc(box / dx))
    grid3d = np.vstack(np.meshgrid(grid1d, grid1d, grid1d)).reshape(3, -1).T

    # get the largest spherical void
    dmax = dx
    for grid_point in grid3d:
        # get the distance (considering minimum image) to the closest atom
        dmin = box
        for atom in zip(frame.x, frame.y, frame.z):
            d = np.linalg.norm([x - box * np.rint(x / box) for x in atom - grid_point])
            dmin = d if d < dmin else dmin

        # get the largest min distance and the positions in the grid
        if dmin > dmax:
            dmax = dmin
            largest_pos = grid_point

    # add a Li atom to the frame
    frame.natoms += 1
    frame.types = np.append(frame.types, "Li")
    frame.x = np.append(frame.x, largest_pos[0])
    frame.y = np.append(frame.y, largest_pos[1])
    frame.z = np.append(frame.y, largest_pos[2])

    # increase the volume and scale all the coordinates
    expand_factor = np.cbrt(box**3 + 16.05) / box
    frame.box *= expand_factor
    frame.x *= expand_factor
    frame.y *= expand_factor
    frame.z *= expand_factor

    return dmax, frame


if __name__ == "__main__":
    box = 10.937456

    # get the last frame of a trajectory
    frames = exma.read_xyz("a-Si64.xyz")
    frame = frames[-1]

    # wrap the frame inside the box
    frame.box = np.full(3, box)
    frame._wrap()

    dmax, frame = lithiation_step(frame, box, 1)

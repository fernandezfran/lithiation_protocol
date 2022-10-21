#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two first steps of the lithiation of amorphous silicon protocol."""
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
        frame_point = exma.core.AtomicSystem(
            natoms=1,
            box=np.full(3, box),
            types=np.array(["Li"]),
            x=np.array([grid_point[0]]),
            y=np.array([grid_point[1]]),
            z=np.array([grid_point[2]]),
        )
        dmin = np.min(exma.distances.pbc_distances(frame_point, frame))

        # get the largest min distance and the positions in the grid
        if dmin > dmax:
            dmax = dmin
            largest_pos = grid_point

    # add a Li atom to the frame
    frame.natoms += 1
    frame.types = np.append(frame.types, "Li")
    frame.x = np.append(frame.x, largest_pos[0])
    frame.y = np.append(frame.y, largest_pos[1])
    frame.z = np.append(frame.z, largest_pos[2])

    # increase the volume and scale all the coordinates
    expand_factor = np.cbrt(box ** 3 + 16.05) / box
    frame.box *= expand_factor
    frame.x *= expand_factor
    frame.y *= expand_factor
    frame.z *= expand_factor

    return dmax, frame

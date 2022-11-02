#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two first steps of the lithiation of amorphous silicon protocol."""
import exma
import numpy as np
import scipy.spatial


def lithiation_step(frame, box):
    """A single step of the lithiation protocol.

    Add a Li atom at the center of the largest spherical void in a structure
    and expand the volume.

    To find the largest spherical void, the centers of the Delaunay
    triangulation are found, which correspond to the vertices of a Voronoi
    diagram, the distance from these points to all the atoms is calculated,
    the smallest one is selected and then the largest of these corresponds to
    the empty sphere with the largest radius.

    Parameters
    ----------
    frame : exma.core.AtomicSystem
        frame of the structure to lithiate and expand

    box : float
        the box size

    Returns
    -------
    tuple
        with a float corresponding with the radius of the sphere and a frame
        (exma.core.AtomicSystem) with the new structure
    """
    # find the voronoi vertices inside the box with pbc
    replicated_frame = exma.io.positions.replicate(frame, [3, 3, 3])
    x = replicated_frame.x
    y = replicated_frame.y
    z = replicated_frame.z

    voronoi = scipy.spatial.Voronoi(np.array((x, y, z)).T)

    vcenter = np.full(3, box)
    mask = (voronoi.vertices >= box) & (voronoi.vertices < 2 * box)
    vertices = [v - vcenter for v, m in zip(voronoi.vertices, mask) if m.all()]

    # get the largest spherical void
    dmax = 1e-6
    for vpoint in vertices:
        # get the distance (considering minimum image) to the closest atom
        frame_point = exma.core.AtomicSystem(
            natoms=1,
            box=np.full(3, box),
            types=np.array(["Li"]),
            x=np.array([vpoint[0]]),
            y=np.array([vpoint[1]]),
            z=np.array([vpoint[2]]),
        )
        dmin = np.min(exma.distances.pbc_distances(frame_point, frame))

        # get the largest min distance and the positions in the grid
        if dmin > dmax:
            dmax = dmin
            largest_pos = vpoint

    # add a Li atom to the frame
    frame.natoms += 1
    frame.types = np.append(frame.types, "Li")
    frame.x = np.append(frame.x, largest_pos[0])
    frame.y = np.append(frame.y, largest_pos[1])
    frame.z = np.append(frame.z, largest_pos[2])

    # increase the volume and scale all the coordinates
    expand_factor = np.cbrt(box**3 + 16.05) / box
    frame.box *= expand_factor
    frame.x *= expand_factor
    frame.y *= expand_factor
    frame.z *= expand_factor

    return dmax, frame

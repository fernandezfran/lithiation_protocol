#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess

import MDAnalysis as mda
from MDAnalysis.transformations import wrap

import numpy as np

from scipy.spatial import Voronoi


class LithiationProtocol:
    def __init__(
        self,
        structure=None,
        nsteps=1,
        nsi=64,
        expansion_factor=16.05,
        xfull=3.75,
        boxsize=10.937456,
    ):
        self.structure = structure
        self.nsteps = nsteps
        self.nsi = nsi
        self.expansion_factor = expansion_factor
        self.xfull = xfull
        self.boxsize = boxsize

    def _get_min_press_frame(self):
        box, press = [], []

        with open("md.out", "r") as md_info:
            line = md_info.readline()

            while line not in (None, "\r", "\n", ""):
                line = line.strip()
                if line.startswith("Pressure"):
                    press.append(line.split()[3])
                elif line.startswith("Lattice vectors"):
                    line = md_info.readline()
                    box.append(line.split()[0])
                line = md_info.readline()

            box = np.array(box, dtype=float)
            press = np.array(press, dtype=float)

        min_press_frame = np.argmin(np.abs(press))

        u = mda.Universe("LixSi64.xyz")
        u.trajectory.dimensions = np.array(3 * [box[min_press_frame]] + 3 * [90.0])
        u.trajectory.add_transformations(wrap(u.atoms))

        return u.trajectory[min_press_frame]

    def _write_gen_format(self, frame):
        with open("LixSi64.gen", "w") as gen:
            gen.write(f" {frame.n_atoms} S\n")
            gen.write(" Si Li\n")
            for i in range(frame.n_atoms):
                atom_type = 1 if i < 64 else 2
                x, y, z = frame.positions[i]
                gen.write(f"{i + 1} {atom_type} {x:.6e} {y:.6e} {z:.6e}\n")
            gen.write("0.0 0.0 0.0\n")
            gen.write(f"{frame.dimensions[0]} 0.0 0.0\n")
            gen.write(f"0.0 {frame.dimensions[1]} 0.0\n")
            gen.write(f"0.0 0.0 {frame.dimensions[2]}\n")

    def _lithiation_step(self, frame, box):
        # scipy does not allow pbc so the system is replicated in all directions
        replicated_frame = exma.io.positions.replicate(frame, [3, 3, 3])
        x = replicated_frame.x
        y = replicated_frame.y
        z = replicated_frame.z

        voronoi = Voronoi(np.array((x, y, z)).T)

        # get the vertices of the voronoi diagram in the central box
        vcenter = np.full(3, box)
        mask = (voronoi.vertices >= box) & (voronoi.vertices < 2 * box)
        vertices = [v - vcenter for v, m in zip(voronoi.vertices, mask) if m.all()]

        # get the largest spherical void
        dmins = []
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
            dmins.append(np.min(exma.distances.pbc_distances(frame_point, frame)))

            del frame_point

        idx = np.argmax(dmins)
        dmax = dmins[idx]
        largest_pos = vertices[idx]

        # add a Li atom to the frame
        frame.natoms += 1
        frame.types = np.append(frame.types, "Li")
        frame.x = np.append(frame.x, largest_pos[0])
        frame.y = np.append(frame.y, largest_pos[1])
        frame.z = np.append(frame.z, largest_pos[2])

        # increase the volume and scale all the coordinates
        expand_factor = np.cbrt(box ** 3 + self.expansion_factor) / box
        frame.box *= expand_factor
        frame.x *= expand_factor
        frame.y *= expand_factor
        frame.z *= expand_factor

        del replicated_frame, x, y, z, voronoi

        return dmax, frame

    def _lithiate_nsteps(self, frame):
        for _ in range(self.nsteps):
            dmax, frame = self._lithiation_step(frame, frame.dimensions[0])
            self.nli_ += 1

        return frame

    def run(self):
        if self.structure is None:
            u = mda.Universe("init/a-Si64.xyz")
            u.trajectory.dimensions = np.array(3 * [self.boxsize] + 3 * [90.0])
            u.trajectory.add_transformations(wrap(u.atoms))
            frame = u.trajectory[-1]

            self.nli_ = 0
            frame = self._lithiate_nsteps(frame)
            self._write_gen_format(frame)
            subprocess.run(["bash", "run.sh"])

        else:
            os.system(f"cp npt/md.{self.structure}.out md.out")
            os.system(f"cp npt/{self.structure}.xyz LixSi64.xyz")

            u = mda.Universe("LixSi64.xyz")
            self.nli_ = np.count_nonzero(u.atoms.types == "LI")

        x = self.nli_ / self.nsi
        while x < self.xfull:
            frame = self._get_min_press_frame()

            os.system(f"mv md.out npt/md.Li{self.nli_}Si64.out")
            os.system(f"mv LixSi64.xyz npt/Li{self.nli_}Si64.xyz")

            frame = self._lithiate_nsteps(frame)

            self._write_gen_format(frame)

            subprocess.run(["bash", "run.sh"])

            x = self.nli_ / self.nsi

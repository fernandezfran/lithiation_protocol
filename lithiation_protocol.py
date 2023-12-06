#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import itertools as it
import os
import subprocess

import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.coordinates.timestep import Timestep
from MDAnalysis.transformations import wrap

import numpy as np

from scipy.spatial import Voronoi


class LithiationProtocol:
    def __init__(
        self,
        restart_from=None,
        nsteps=1,
        nsi=64,
        expansion_factor=16.05,
        xfull=3.75,
        box_size=10.937456,
    ):
        self.restart_from = restart_from
        self.nsteps = nsteps
        self.nsi = nsi
        self.expansion_factor = expansion_factor
        self.xfull = xfull
        self.box_size = box_size

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

    def _lithiation_step(self, frame):
        box = frame.dimensions[:3]
        images = box * list(it.product((-1, 0, 1), repeat=3))
        positions = np.concatenate([frame.positions + image for image in images])

        voronoi = Voronoi(positions)

        mask = (voronoi.vertices >= 0.0) & (voronoi.vertices < box)
        vertices = [v for v, m in zip(voronoi.vertices, mask) if m.all()]

        distances = distance_array(positions, vertices, box=box, backend="OpenMP")
        dmins = [np.min(dist[dist > 0] for dist in distances)]

        sf = np.cbrt(box**3 + self.expansion_factor) / box

        ts = Timestep(self.nsi + self.nli_)
        ts.dimensions = np.concatenate((sf * box, np.full(3, 90)))
        new_li = [vertices[np.argmax(dmins)]]
        ts.positions = sf * np.concatenate((frame.positions, new_li))

        return ts

    def _lithiate_nsteps(self, frame):
        for _ in range(self.nsteps):
            self.nli_ += 1
            frame = self._lithiation_step(frame)

        return frame

    def run(self):
        if self.restart_from is None:
            u = mda.Universe("init/a-Si64.xyz")
            u.trajectory.dimensions = np.array(3 * [self.box_size] + 3 * [90.0])
            u.trajectory.add_transformations(wrap(u.atoms))
            frame = u.trajectory[-1]

            self.nli_ = 0
            frame = self._lithiate_nsteps(frame)
            self._write_gen_format(frame)
            subprocess.run(["bash", "run.sh"])

        else:
            os.system(f"cp npt/md.{self.restart_from}.out md.out")
            os.system(f"cp npt/{self.restart_from}.xyz LixSi64.xyz")

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

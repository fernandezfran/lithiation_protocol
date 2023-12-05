#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Full lithiation of amorphous silicon using DFTB+.

We start with an amorphized silicon structure and follow the next protocol:

1. add a Li atom at the center of the largest spherical void,
2. increase the volume and scale the coordinates,
3. perform a local LBFGS minimization,
4. simulate an NPT molecular dynamics,
5. select the frame with minimum absolute pressure,
6. with x defined as number of Li atoms per Si atoms if x < 3.75 goto point 1
else finish.

To make the process faster, till 16 atoms of lithium can be added at a time
and expanding in each one of them, this was checked with DFT which does not
alter the results.
"""
import os
import subprocess

import MDAnalysis as mda

import numpy as np

from lithiation_step import lithiation_step


class LithiationProtocol:
    """Lithiation protocol.

    Parameters
    ----------
    structure : str, default=None
        name of the structure to restart from

    nsteps : int, default=1
        number of simultaneous lithium insertions
    """

    def __init__(self, structure=None, nsteps=1, nsi=64):
        self.structure = structure
        self.nsteps = nsteps
        self.nsi = nsi

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
        frame = u.trajectory[min_press_frame]
        frame.dimensions = np.array(3 * [box[min_press_frame]] + 3 * [90.0])

        # falta hacer un wrap de las posiciones acá

        return frame

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


    def _lithiate_nsteps(self, frame):
        # add nsteps Li atoms and expand frame
        for _ in range(self.nsteps):
            dmax, frame = lithiation_step(frame, frame.dimensions[0])
            self.nli += 1

        return frame

    def run(self):
        """Run the full lithiation or restart from an structure."""
        if self.structure is not None:
            # get the structure and restart the lithiation
            os.system(f"cp npt/md.{self.structure}.out md.out")
            os.system(f"cp npt/{self.structure}.xyz LixSi64.xyz")

            # initial number of Li atoms
            frames = exma.read_xyz(f"LixSi64.xyz")
            frame = frames[-1]
            self.nli = frame._natoms_type(frame._mask_type("Li"))

        else:
            # initializate lithiation
            u = mda.Universe("init/a-Si64.xyz")
            frame = u.trajectory[-1]
            frame.dimensions = np.array(3 * [10.937456] + 3 * [90.0])

            # falta hacer un wrap de las posiciones acá

            # add nsteps Li atom and expand frame
            self.nli = 0
            frame = self._lithiate_nsteps(frame)
            self._write_gen_format(frame)
            subprocess.run(["bash", "run.sh"])

        x = self.nli / self.nsi
        while x < 3.75:
            frame = self._get_min_press_frame()

            # save files of interest
            os.system(f"mv md.out npt/md.Li{self.nli}Si64.out")
            os.system(f"mv LixSi64.xyz npt/Li{self.nli}Si64.xyz")

            frame = self._lithiate_nsteps(frame)

            # write genFormat file for dftb+
            self._write_gen_format(frame)

            # run LBFGS minimization and Berendsen NPT
            subprocess.run(["bash", "run.sh"])

            x = self.nli / self.nsi

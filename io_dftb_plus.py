#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd


def write_gen_format(frame, filename):
    """Write a .gen format file with a xyz frame.

    Parameters
    ----------
    frame : exma.core.AtomicSystem
        frame of the structure to write as genFormat

    filename : str
        name of the .gen file
    """
    with open(filename, "w") as gen:
        gen.write(f" {frame.natoms} S\n")
        gen.write(" Si Li\n")
        for i, (t, x, y, z) in enumerate(zip(frame.types, frame.x, frame.y, frame.z)):
            atom_type = 1 if t == "Si" else 2
            gen.write(f"{i + 1} {atom_type} {x:.6e} {y:.6e} {z:.6e}\n")
        gen.write("0.0 0.0 0.0\n")
        gen.write(f"{frame.box[0]} 0.0 0.0\n")
        gen.write(f"0.0 {frame.box[1]} 0.0\n")
        gen.write(f"0.0 0.0 {frame.box[2]}\n")

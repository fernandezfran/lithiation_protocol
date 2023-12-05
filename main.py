#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Lithiation protocol of Si."""

import argparse

import lithiation_protocol


def main():
    parser = argparse.ArgumentParser(
        description="Lithiate an amorphous structure, by default from the "
        "beginning but can also be restarted from a given structure."
    )

    for flag, description, default in zip(
        ["restart-from", "nsteps", "natoms", "expansion-factor"],
        [
            "restart from a given structure RESTART_FROM, e.g. Li55Si64",
            "number of simultaneous lithium insertions, e.g. 3",
            "number of atoms in the initial amorphous structure, e.g. 64",
            "the volume expansion of adding a lithium atom in the structure",
        ],
        [None, 1, 16.05, 64],
    ):
        parser.add_argument(f"--{flag}", help=description, default=default)

    args = parser.parse_args()

    lithiation_protocol.LithiationProtocol(**args.__dict__).run()


if __name__ == "__main__":
    main()

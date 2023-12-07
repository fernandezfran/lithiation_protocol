#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import lithiation_protocol


def main():
    parser = argparse.ArgumentParser(
        description="Lithiation protocol for an amorphous Si structure"
    )

    for flag, description, default, dtype in zip(
        [
            "restart-from",
            "nsteps",
            "expansion-factor",
            "nsi",
            "xfull",
            "box-size",
        ],
        [
            "restart from a given structure RESTART_FROM, e.g. Li55Si64",
            "number of simultaneous lithium insertions, e.g. 3",
            "the volume expansion of adding a lithium atom in the structure",
            "number of atoms of Si in the initial amorphous structure, e.g. 64",
            "the maximum x value",
            "the initial size of the box",
        ],
        [None, 1, 16.05, 64, 3.75, 10.937456],
        [str, int, float, int, float, float],
    ):
        parser.add_argument(f"--{flag}", help=description, default=default, type=dtype)

    args = parser.parse_args()

    lp = lithiation_protocol.LithiationProtocol(**vars(args))
    lp.run()


if __name__ == "__main__":
    main()

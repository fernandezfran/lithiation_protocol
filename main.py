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

    parser.add_argument(
        "--restart-from",
        help="restart from a given structure RESTART_FROM, e.g. Li55Si64",
        default=None,
    )

    parser.add_argument(
        "--nsteps",
        help="number of simultaneous lithium insertions, e.g. 3",
        default=1,
    )

    args = parser.parse_args()

    lithiation_protocol.LithiationProtocol(**args.__dict__).run()


if __name__ == "__main__":
    main()

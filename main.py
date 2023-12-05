#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Lithiation protocol of Si."""

import argparse

import lithiation


def main():
    parser = argparse.ArgumentParser(
        description="Lithiate an amorphous structure, "
        "by default from the beginning but can also be restarted from a given "
        "structure. You have also different options to analyze the lithiation "
        " once it was performed."
    )

    parser.add_argument(
        "--restart-from",
        help="restart from a given structure RESTART_FROM, e.g. Li55Si64",
    )

    args = parser.parse_args()

    full = True
    for v in args.__dict__.values():
        if v:
            full = False

    if full or args.restart_from:
        lithiation.Lithiation(structure=args.restart_from, nsteps=1).run()


if __name__ == "__main__":
    main()

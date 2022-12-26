#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Different mechanical statistical analyses of the lithiation."""
import os

import exma

import matplotlib.pyplot as plt
import matplotlib.colors

import numpy as np

import pandas as pd

from io_dftb_plus import read_md_out


NSI = 64
NLIMAX = 3.75 * NSI


def _color_fader(mix, first_color="#1f77b4", second_color="#2ca02c"):
    c1 = np.array(matplotlib.colors.to_rgb(first_color))
    c2 = np.array(matplotlib.colors.to_rgb(second_color))
    return matplotlib.colors.to_hex((1 - mix) * c1 + mix * c2)


def _colormap():
    colors = [_color_fader(v) for v in np.linspace(0, 1, num=100)]
    return matplotlib.colors.ListedColormap(colors)


def _read_trajectories(path="npt/"):
    """Read all the xyz files and md outs in `npt` directory."""
    trajectories = [
        exma.read_xyz(path + fname)
        for fname in os.listdir(path)
        if fname.endswith(".xyz")
    ]

    natoms = np.array([traj[0].natoms for traj in trajectories])
    args = np.argsort(natoms)
    trajectories = np.array(trajectories)[args]

    info = []
    for traj in trajectories:
        nli = traj[0]._natoms_type(traj[0]._mask_type("Li"))
        info.append(read_md_out(filename=path + f"md.Li{nli}Si{NSI}.out"))

    return trajectories, info


def rdf_plot(each=5, central="Si", interact="Si", save=False):
    """Plot the RDF of Si-Si for each lithium concentration."""
    trajectories, info = _read_trajectories()

    fig, ax = plt.subplots()
    cmap = _colormap()

    for k, (frames, thermo) in enumerate(zip(trajectories, info)):

        if k % each == 0:
            nli = frames[0]._natoms_type(frames[0]._mask_type("Li"))

            for box, frame in zip(thermo["xbox"].values, frames):
                frame.box = np.full(3, box)

            rdf = exma.rdf.RadialDistributionFunction(
                frames, type_c=central, type_i=interact, rmax=5.0
            )
            rdf.calculate()
            rdf.plot(ax=ax, plot_kws={"color": cmap(nli / NLIMAX)})

    ax.set_xlabel(r"r [$\AA$]")
    ax.set_ylabel(f"RDF {central}-{interact}")

    ax.set_xlim((1.5, 5))

    ax.grid(linestyle=":")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    clb = fig.colorbar(sm)
    clb.ax.set_ylabel("Percent lithium concentration")

    fig.tight_layout()
    if save:
        fig.savefig(f"res/rdf_{central}-{interact}.png", dpi=600)
    plt.show()


def fvc_plot(save=False):
    """Plot the fractional volume change during the lithiation."""
    trajectories, info = _read_trajectories()

    x_values, volume, err_volume = [], [], []
    for frames, thermo in zip(trajectories, info):
        nli = frames[0]._natoms_type(frames[0]._mask_type("Li"))

        vx = [box ** 3 for box in thermo["xbox"].values]

        x_values.append(nli / NSI)
        volume.append(np.mean(vx))
        err_volume.append(np.std(vx))

    df = pd.DataFrame(
        {
            "x": x_values,
            "natoms_a": np.full(len(volume), NSI),
            "volume": volume,
            "err_volume": err_volume,
        }
    )
    df = exma.electrochemistry.fractional_volume_change(df, NSI, 10.937456 ** 3)

    fig, ax = plt.subplots()

    ax.errorbar(df["x"], df["fvc"], yerr=df["err_fvc"], marker="o", linestyle=":")

    fvc_chevrier = lambda x: 0.786647 * x + 0.001547
    xeval = np.array([0, 3.75])
    ax.plot(xeval, fvc_chevrier(xeval))

    ax.set_xlabel(r"$x$ in Li$_x$Si")
    ax.set_ylabel("Fractional volume change")

    ax.set_xlim((0, 3.75))

    ax.grid(linestyle=":")

    fig.tight_layout()
    if save:
        fig.savefig("res/fvc.png", dpi=600)
    plt.show()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Different mechanical statistical analyses of the lithiation."""
import exma
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import pandas as pd

from io_dftb_plus import read_md_out


NSI = 64

NORDBLUE = "#5E81AC"
NORDGREEN = "#A3BE8C"


def color_fader(mix):
    c1 = np.array(matplotlib.colors.to_rgb(NORDBLUE))
    c2 = np.array(matplotlib.colors.to_rgb(NORDGREEN))
    return matplotlib.colors.to_hex((1 - mix) * c1 + mix * c2)


def colormap():
    colors = [color_fader(v) for v in np.linspace(0, 1, num=100)]
    return matplotlib.colors.ListedColormap(colors)


def rdf_plot(nlis, path="npt/"):
    """Plot the RDF of Si-Si for each lithium concentration."""
    fig, ax = plt.subplots()
    cmap = colormap()

    for nli in nlis:
        prefix = f"Li{nli}Si{NSI}"

        thermo = read_md_out(path + "md." + prefix + ".out")
        frames = exma.read_xyz(path + prefix + ".xyz")

        for box, frame in zip(thermo["xbox"].values, frames):
            frame.box = np.full(3, box)

        rdf = exma.rdf.RadialDistributionFunction(
            frames, type_c="Si", type_i="Si", rmax=5.0
        )
        rdf.calculate()
        rdf.plot(ax=ax, plot_kws={"color": cmap(nli / (3 + max(nlis)))})

    ax.set_xlabel(r"r [$\AA$]")
    ax.set_ylabel("RDF Si-Si")

    ax.set_xlim((1.5, 5))

    ax.grid(linestyle=":")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    clb = fig.colorbar(sm)
    clb.ax.set_ylabel("Percent lithium concentration")

    fig.tight_layout()
    # fig.savefig(f"res/rdf.png", dpi=600)
    plt.show()


def fvc_plot(nlis, path="npt/"):
    """Plot the fractional volume change during the lithiation."""
    x_values = np.asarray(nlis) / NSI
    natoms_a = np.full(len(nlis), NSI)

    volume, err_volume = [], []
    for x, nli in zip(x_values, nlis):
        prefix = f"Li{nli}Si{NSI}"

        thermo = read_md_out(path + "md." + prefix + ".out")

        vx = [box**3 for box in thermo["xbox"].values]

        volume.append(np.mean(vx))
        err_volume.append(np.std(vx))

    df = pd.DataFrame(
        {"x": x_values, "natoms_a": natoms_a, "volume": volume, "err_volume": volume}
    )
    df = exma.electrochemistry.fractional_volume_change(df, NSI, 10.937456**3)

    fig, ax = plt.subplots()
    ax.errorbar(df["x"], df["fvc"], yerr=df["err_fvc"])

    ax.set_xlabel(r"$x$ in Li$_x$Si")
    ax.set_ylabel("Fractional volume change")

    ax.grid(linestyle=":")

    fig.tight_layout()
    #fig.savefig(f"res/fvc.png", dpi=600)
    plt.show()


if __name__ == "__main__":
    nlis = np.arange(1, 26, 3)
    rdf_plot(nlis)
    fvc_plot(nlis)

# Lithiation protocol of Si

[![PRB](https://img.shields.io/badge/PhysRevB-108.144201-b31033)](https://doi.org/10.1103/PhysRevB.108.144201)
[![MIT](https://img.shields.io/badge/License-MIT-ffde57)](https://github.com/fernandezfran/lithiation_protocol/blob/main/LICENSE)

We start with an amorphous silicon (a-Si) structure and follow the next protocol:

1. add a Li atom at the center of the largest spherical void,
    to find the largest spherical void, the centers of the Delaunay
    triangulation are found, which correspond to the vertices of a Voronoi
    diagram, the distance from these points to all the atoms is calculated,
    the smallest one is selected and then the largest of these corresponds to
    the empty sphere with the largest radius.

2. increase the volume and scale the coordinates,
    this is done to follow the experimental expansion of the system.

3. perform a local LBFGS minimization,
    with DFTB+ software.

4. simulate an NPT molecular dynamics,
    with DFTB+ software.

5. select the frame with minimum absolute pressure,

6. with x defined as number of Li atoms per Si atoms if x is less than 3.75 goto point 1 else finish.

This protocol is slightly similar to the one proposed by 
[Chevrier and Dahn](https://doi.org/10.1149/1.3111037), also for the lithiation 
of a-Si. To make the process faster, multiple atoms of lithium can be added at 
a time and expanding the volume and the coordinates in each one of them. The 
authors checked with DFT that a step of x=0.25 in Li<sub>x</sub>Si does not alter 
the results, which corresponds with 16 Li atoms for the a-Si64 initial structure.


## Requirements

Python3.8+ and [MDAnalysis](https://www.mdanalysis.org/) (which also install 
NumPy and SciPy):

```
pip install -U MDAnalysis
```

You also need a [DFTB+](https://github.com/dftbplus/dftbplus) executable and the 
set of parameters for the LiSi interaction in the `params` directory (these 
parameters can be downloaded 
[here](https://github.com/alexispaz/DFTB_LiSi/tree/main/lisi)). Also, `hsd` input
files for DFTB+ LBFGS minimization and NPT simulation must be provided in the 
respective directories.


## Usage

```
$ python3 main.py --help
usage: main.py [-h] [--restart-from RESTART_FROM] [--nsteps NSTEPS] [--expansion-factor EXPANSION_FACTOR] [--nsi NSI] [--xfull XFULL]
               [--box-size BOX_SIZE]

Lithiate an amorphous structure, by default from the beginning but can also be restarted from a given structure.

optional arguments:
  -h, --help            show this help message and exit
  --restart-from RESTART_FROM
                        restart from a given structure RESTART_FROM, e.g. Li55Si64
  --nsteps NSTEPS       number of simultaneous lithium insertions, e.g. 3
  --expansion-factor EXPANSION_FACTOR
                        the volume expansion of adding a lithium atom in the structure
  --nsi NSI             number of atoms of Si in the initial amorphous structure, e.g. 64
  --xfull XFULL         the maximum x value
  --box-size BOX_SIZE   the initial size of the box
```

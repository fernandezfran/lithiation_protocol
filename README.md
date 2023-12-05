# Lithiation protocol of Si

We start with an amorphous silicon (a-Si) structure and follow the next protocol:

1. add a Li atom at the center of the largest spherical void,
2. increase the volume and scale the coordinates,
3. perform a local LBFGS minimization,
4. simulate an NPT molecular dynamics,
5. select the frame with minimum absolute pressure,
6. with x defined as number of Li atoms per Si atoms if x is less than 3.75 goto point 1 else finish.

This protocol is slightly similar to the one proposed by 
[Chevrier and Dahn](https://doi.org/10.1149/1.3111037), also for the lithiation 
of a-Si. To make the process faster, multiple atoms of lithium can be added at 
a time and expanding in each one of them. The authors checked with DFT that a 
step of x=0.25 in Li<sub>x</sub>Si does not alter the results, which would
correspond to 16 Li atoms for the a-Si64 structure.


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
usage: main.py [-h] [--restart-from RESTART_FROM] [--nsteps NSTEPS] [--natoms NATOMS] [--expansion-factor EXPANSION_FACTOR]

Lithiate an amorphous structure, by default from the beginning but can also be restarted from a given structure.

optional arguments:
  -h, --help            show this help message and exit
  --restart-from RESTART_FROM
                        restart from a given structure RESTART_FROM, e.g. Li55Si64
  --nsteps NSTEPS       number of simultaneous lithium insertions, e.g. 3
  --natoms NATOMS       number of atoms in the initial amorphous structure
  --expansion-factor EXPANSION_FACTOR
                        the volume expansion of adding a lithium atom in the structure

```

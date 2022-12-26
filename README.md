# Lithiation of amorphous silicon

We start with an amorphized silicon structure and follow the next protocol:

1. add a Li atom at the center of the largest spherical void,
2. increase the volume and scale the coordinates,
3. perform a local LBFGS minimization,
4. simulate 10ps with Berendsen NPT,
5. select the frame with minimum absolute pressure,
6. with x defined as number of Li atoms per Si atoms if x is less than 3.75 goto point 1 else finish.

To make the process faster, multiple atoms of lithium are added at a time and 
expanding in each one of them. It was checked with DFT that a step of x=0.25 in
Li$_x$Si does not alter the results, which would correspond to 16 atoms here.

This protocol is slightly similar to the one proposed by 
[Chevrier and Dahn](https://doi.org/10.1149/1.3111037), also for the lithiation 
of amorphous silicon.


## Requirements

Python3.8+ and the required libraries:

```bash
pip install -r requirements.txt
```

You also need a [DFTB+](https://github.com/dftbplus/dftbplus) executable. In case
the minimizations and equilibrations want to be done with other software you can 
modify the `run.sh` file to do these two steps in another way.


## Usage

```bash
$ python3 main.py --help
usage: main.py [-h] [--restart-from RESTART_FROM] [--fvc] [--rdf [RDF]] [--central [CENTRAL]] [--interact [INTERACT]] [-s]

Lithiate an amorphous structure, by default from the beginning but can also be restarted from a given structure. You have also different
options to analyze the lithiation once it was performed.

optional arguments:
  -h, --help            show this help message and exit
  --restart-from RESTART_FROM
                        restart from a given structure RESTART_FROM, e.g. Li55Si64
  --fvc                 fractional volume change calculation
  --rdf [RDF]           Si-Si rdf specifying how many structures to skip for the plot, by default 5
  --central [CENTRAL]   central atom type for rdf calculation, by default Si
  --interact [INTERACT]
                        interact atom type for rdf calculation, by default Si
  -s, --save            true to save the png figures if created
```

For example, if you want to start the lithiation from scratch:
```bash
$ python3 main.py 
```
But if you want to restart from the structure, e.g. Li17Si64, then
```bash
$ python3 main.py --restart-from Li17Si64
```

Once the lithiation has been run, it can be analyzed structurally, for example, 
the fractional volumetric change 
```bash
$ python3 main.py --fvc
```
or the rdf
```bash
$ python3 main.py --rdf
```
if you want to change the default Si-Si interaction you can use the following flags
```bash
$ python3 main.py --rdf --central Si --interact Li
```
if you want to save the plots you can use the flag `-s` at the end.

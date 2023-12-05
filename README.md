# Lithiation protocol of Si

We start with an amorphized silicon structure and follow the next protocol:

1. add a Li atom at the center of the largest spherical void,
2. increase the volume and scale the coordinates,
3. perform a local LBFGS minimization,
4. simulate an NPT molecular dynamics,
5. select the frame with minimum absolute pressure,
6. with x defined as number of Li atoms per Si atoms if x is less than 3.75 goto point 1 else finish.

This protocol is slightly similar to the one proposed by 
[Chevrier and Dahn](https://doi.org/10.1149/1.3111037), also for the lithiation 
of amorphous silicon. To make the process faster, multiple atoms of lithium 
are added at a time and expanding in each one of them. The authors checked 
with DFT that a step of x=0.25 in Li$_x$Si does not alter the results, 
which would correspond to 16 atoms here.


## Requirements

Python3.8+ and the required libraries:

```
pip install -r requirements.txt
```

You also need a [DFTB+](https://github.com/dftbplus/dftbplus) executable. 

The set of parameters for the LiSi interaction can be downloaded 
[here](https://github.com/alexispaz/DFTB_LiSi/tree/main/lisi) and you must copy 
them to the `params` directory.

In case the minimizations and equilibrations want to be done with other 
software and/or potential you can modify the `run.sh` file to do these two 
steps in another way.


## Usage

```
$ python3 main.py --help
usage: main.py [-h] [--restart-from RESTART_FROM]

Lithiate an amorphous structure, by default from the beginning but can also be restarted from a given structure. 

optional arguments:
  -h, --help            show this help message and exit
  --restart-from RESTART_FROM
                        restart from a given structure RESTART_FROM, e.g. Li55Si64
```

For example, if you want to start the lithiation from scratch:
```
$ python3 main.py 
```
But if you want to restart from the structure, e.g. Li17Si64, then
```
$ python3 main.py --restart-from Li17Si64
```

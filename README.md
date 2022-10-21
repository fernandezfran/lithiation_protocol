# Lithiation of amorphous silicon

Slightly similar [Chevrier and Dahn protocol](https://doi.org/10.1149/1.3111037)
for lithiation of a-Si.

We start with an amorphous Si structure and follow the protocol:
1. add a Li atom at the center of the largest spherical void,
2. increase the volume and scale the coordinates,
3. perform a local LBFGS minimization
4. simulate 10ps with Berendsen NPT,
5. select the frame with minimum absolute pressure,
6. with x defined as number of Li atoms per Si atoms if x is less than 3.75 goto point 1 else finish.

To make the process faster, 3 atoms of lithium are added at a time and expanding
in each one of them, this was checked with DFT which does not alter the results.

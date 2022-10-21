# Lithiation of amorphous silicon

Slightly similar [Chevrier and Dahn protocol](https://doi.org/10.1149/1.3111037)
for lithiation of a-Si.

We start with an amorphous Si structure and follow the protocol:
1. add a Li atom at the center of the largest spherical void,
2. increase the volume and scale the coordinates,
3. simulate 5ps with Berendsen NPT,
4. select the frame with minimum absolute pressure,
5. goto point 1 if x is less than 3.75 else finish.

To make the process faster, 3 atoms of lithium are added at a time and expanding
in each one of them, this was checked with DFT which does not alter the results.

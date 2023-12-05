#!bin/bash

# move the gen file to the lbfgs dir
mv LixSi64.gen lbfgs/

# goto lbfgs dir
cd lbfgs/

# run lbfgs minimization
dftb+

# copy the minimized structure to the npt dir 
cp LixSi64.gen ../npt/

# goto npt dir
cd ../npt/

# run npt molecular dynamics simulation
dftb+

# copy the generated files to the main dir
cp md.out ../
cp LixSi64.xyz ../

# goback main dir
cd ../

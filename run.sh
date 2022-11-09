#!bin/bash
mv LixSi64.gen lbfgs/
cd lbfgs/
dftb+
cp LixSi64.gen ../npt/
cd ../npt/
dftb+
cp md.out ../
cp LixSi64.xyz ../
cd ../

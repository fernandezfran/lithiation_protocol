#!bin/bash
mv LixSi64.gen dftb+min/
cd dftb+min/
dftb+
cp LixSi64.gen ../dftb+npt/
cd ../dftb+npt/
dftb+
cp md.out ../
cp LixSi64.xyz ../
cd ../

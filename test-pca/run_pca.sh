#!/bin/bash

suffix=printcov
cp ../pca_$suffix/pca_$suffix .
./pca_$suffix -irf ref.txt -itj trj.txt -imw mass.txt -ipr pca.in -ors rmsd.dat -oav ave.dat -ova val.dat -ocf cfl.dat -ove vec.dat -cov cov.dat -ocd covd.dat

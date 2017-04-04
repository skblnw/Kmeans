#!/bin/bash

suffix=newinitial
ifort -p -o Kmeans_$suffix ../src/RC_Kmeans_$suffix.f90
mkdir -p $suffix
cd $suffix
rm -f *
cp ../barph/cov* .
../Kmeans_$suffix

#!/bin/bash

suffix=newinitial
ifort -p -o Kmeans_$suffix ../src/RC_Kmeans_$suffix.f90
rm -rf $suffix
mkdir -p $suffix
cd $suffix
cp ../barph/cov* .
../Kmeans_$suffix

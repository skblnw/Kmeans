#!/bin/bash

suffix=newinitial
ifort -p -o Kmeans_$suffix ../src/RC_Kmeans/RC_Kmeans_$suffix.f90
if true; then
rm -rf $suffix
mkdir -p $suffix
cd $suffix
cp ../hsp/cov* .
../Kmeans_$suffix
fi

#!/bin/bash

suffix=newinitial
ifort -p -o Kmeans_$suffix ../RC_Kmeans_$suffix.f90
mkdir -p $suffix
cd $suffix
rm -f *
cp ../cov* .
../Kmeans_$suffix

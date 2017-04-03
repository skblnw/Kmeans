#!/bin/bash

suffix=newinitial
ifort -p -o Kmeans_$suffix ../RC_Kmeans_$suffix.f90
cd $suffix
rm -f *.txt
../Kmeans_$suffix

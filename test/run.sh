#!/bin/bash

suffix=newinitial
ifort -p -o Kmeans_$suffix ../RC_Kmeans_$suffix.f90
#./Kmeans_$suffix

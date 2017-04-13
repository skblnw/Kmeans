#!/usr/bin/env python3

import argparse
import numpy as np
import sys
sys.path.append("/home/kevin/scripts/mkpy")
import mkpy

parser = argparse.ArgumentParser()
parser.add_argument("ced", help="Input coved.dat")
parser.add_argument("cod", help="Input covd.dat")
parser.add_argument("natm", type=int, help="Number of atoms")
parser.add_argument("ncg", type=int, help="Desired number of CG sites")

args = parser.parse_args()


print("Reading " + str(args.ced))
with open(args.ced,'r') as f:
    C = []
    for line in f:
        line = line.split()
        if line:
            line = [float(i) for i in line]
            C.extend(line)
C = np.array(C)
C = mkpy.squareform_diagfill(C)

print("Reading " + str(args.cod))
with open(args.cod,'r') as f:
    D = []
    for line in f:
        line = line.split()
        if line:
            line = [float(i) for i in line]
            D.extend(line)
D = np.array(D)
D = mkpy.squareform_diagfill(D)

A = np.zeros(C.shape,dtype=C.dtype)
for ii in np.arange(C.shape[0]):
    for jj in np.arange(C.shape[0]):
        A[ii,jj] = C[ii,ii] - C[ii,jj]*2 + C[jj,jj]

B = np.zeros(D.shape,dtype=D.dtype)
for ii in np.arange(D.shape[0]):
    for jj in np.arange(D.shape[0]):
        B[ii,jj] = D[ii,ii] - D[ii,jj]*2 + D[jj,jj]


list = []
for atm in range(args.natm):
    cor = A[atm,:]
    ind = [atm]
    for cg in range(args.ncg):
        ind_tmp = np.where( cor == cor.max() )[0][0]
        ind.extend([ind_tmp])
        cor = np.amin(A[ind,:], axis=0)
    ind = sorted(ind)
    list.append([ind])
print(list)
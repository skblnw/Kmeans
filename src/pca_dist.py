#!/usr/bin/env python3

import argparse
import numpy as np
import sys
sys.path.append("/home/kevin/scripts/mkpy")
import mkpy

parser = argparse.ArgumentParser()
parser.add_argument("ced", help="Input coved.dat")
parser.add_argument("cod", help="Input covd.dat")

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
print(C.shape)
C = mkpy.squareform_diagfill(C)
print(C.shape)

print("Reading " + str(args.cod))
with open(args.cod,'r') as f:
    D = []
    for line in f:
        line = line.split()
        if line:
            line = [float(i) for i in line]
            D.extend(line)

D = np.array(D)
print(D.shape)
D = mkpy.squareform_diagfill(D)
print(D.shape)

ii = np.arange(C.shape[0])
jj = np.arange(C.shape[0])

A = np.zeros(C.shape,dtype=C.dtype)
for ii in np.arange(C.shape[0]):
    for jj in np.arange(C.shape[0]):
        A[ii,jj] = C[ii,ii] - C[ii,jj]*2 + C[jj,jj]
print(A.shape)

B = np.zeros(D.shape,dtype=D.dtype)
for ii in np.arange(D.shape[0]):
    for jj in np.arange(D.shape[0]):
        B[ii,jj] = D[ii,ii] - D[ii,jj]*2 + D[jj,jj]
print(B.shape)

print(C[0:3,0:3], "\n", A[0:3,0:3])
print(D[0:3,0:3], "\n", B[0:3,0:3])
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
parser.add_argument("lamb", type=float, help="Value of Lamda")

args = parser.parse_args()
lamb = args.lamb

# Read coved.dat
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

# Read covd.dat
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

# Compute A from C
A = np.zeros(C.shape,dtype=C.dtype)
for ii in np.arange(C.shape[0]):
    for jj in np.arange(C.shape[0]):
        A[ii,jj] = C[ii,ii] - C[ii,jj]*2 + C[jj,jj]

# Compute B from D
B = np.zeros(D.shape,dtype=D.dtype)
for ii in np.arange(D.shape[0]):
    for jj in np.arange(D.shape[0]):
        B[ii,jj] = D[ii,ii] - D[ii,jj]*2 + D[jj,jj]

# Compute F from A+B
F = A + lamb * B

# Iniitialization
initial_list = []
for atm in range(args.natm):
    cor = F[atm,:]
    ind = [atm]
    for cg in range(args.ncg-1):
        tmp = np.where( cor == cor.max() )[0]
        ind.extend(tmp)
        cor = np.amin(F[ind,:], axis=0)
    ind = sorted(ind)
    initial_list.append(ind)

# Sort and unique
a = np.array(initial_list)
b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
_, idx = np.unique(b, return_index=True)
unique_initial_list = a[idx]

# Update CG labels using WCSS
# WCSS(C) = Cii - 2/S * Cij + 1/S^2 * Cjj
# All wcss component matrix (wii, wij, wjj) are {natm X ncg} matrix
for index_list in unique_initial_list[:1,:]:
    print(index_list)
    atm_wcss = F[:,index_list]
    print(atm_wcss.shape)
    
    atm_cglabel = np.zeros((args.natm,))
    for ii in range(100):
        # Compare and print atom labeled with CG index
        # Break if converged
        atm_cglabel_new = np.argmin(atm_wcss, axis=1)
        if np.array_equal(atm_cglabel_new, atm_cglabel):
            print("Finished")
            break
        else:
            print("Iteration ", str(ii))
            atm_cglabel = atm_cglabel_new
            print(atm_cglabel)
        atm_Fvalue = np.amin(atm_wcss, axis=1)
        #print(atm_Fvalue)
        
        # Calculate Cjj
        cg_atmlabel = []
        atm_wcss = []
        for cg in range(args.ncg):
            atmlabel = np.where( atm_cglabel == cg)[0]
            cg_atmlabel.append(atmlabel)
            wjj_in_cluster = 0.0
            for j1 in atmlabel:
                for j2 in atmlabel:
                    wjj_in_cluster += C[j1, j2] + lamb * D[j1, j2]
            wjj_in_cluster /= atmlabel.size * atmlabel.size
            atm_wcss.append(wjj_in_cluster)
        atm_wcss = np.tile(np.array(atm_wcss), (args.natm,1))
        #print(atm_wcss[0,:])
        
        # Calculate Cij and add to the atm_wcss matrix
        for cg in range(args.ncg):
            atmlabel = np.where( atm_cglabel == cg)[0]
            wij_in_cluster = 0.0
            for jj in atmlabel:
                wij_in_cluster += C[:, jj] + lamb * D[:, jj]
            wij_in_cluster = wij_in_cluster * 2 / atmlabel.size
            atm_wcss[:,cg] -= wij_in_cluster
        #print(atm_wcss[0,:])
        
        # Calculate Cii and add to the atm_wcss matrix
        wii_in_cluster = C.diagonal() + lamb * D.diagonal()
        wii_in_cluster = np.tile(wii_in_cluster.T, (args.ncg,1)).T
        atm_wcss += wii_in_cluster
        #print(atm_wcss[0,:])
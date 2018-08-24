#!/usr/bin/env python3

import argparse
import numpy as np
import sys
sys.path.append("/home/kevin/scripts/mkpy")
sys.path.append(".")
import mkpy

# FUNCTION DECLARATIONS

def calc_residual_wcss(C, ncg, atm_cglabel):
    # Calculate WCSS per CG sets
    # w(xx)_in_cluster: { (cluster_size) X 1 }
    # css_in_cluster:   { (cluster_size) X 1 }
    cg_wcss = []
    atm_wii = C.diagonal()
    for cg in range(ncg):
        atmlabel = np.where( atm_cglabel == cg)[0]
        wjj_in_cluster = 0.0
        for j1 in atmlabel:
            for j2 in atmlabel:
                wjj_in_cluster += C[j1, j2]
        wjj_in_cluster /= atmlabel.size * atmlabel.size
        css_in_cluster = np.tile(wjj_in_cluster, (atmlabel.size,1))
        
        wij_in_cluster = 0.0
        for jj in atmlabel:
            wij_in_cluster += C[atmlabel, jj]
        wij_in_cluster = wij_in_cluster * 2 / atmlabel.size
        wij_in_cluster = wij_in_cluster.reshape((atmlabel.size,1))
        css_in_cluster -= wij_in_cluster
        
        wii_in_cluster = atm_wii[atmlabel]
        wii_in_cluster = wii_in_cluster.reshape((atmlabel.size,1))
        css_in_cluster += wii_in_cluster
        
        cg_wcss.append(css_in_cluster.sum())
        
    return np.array(cg_wcss)

def calc_wcss_per_atom(C, D, lamb, ncg, atm_cglabel):
    #print("Iteration ", str(ii))
    
    # atm_wcss: { natm X ncg }
    # atm_wjj: { 1    X ncg } -> tile to natm rows
    # atm_wij: { natm X 1   } -> substract from each column of atm_wcss
    # atm_wii: { natm X 1   } -> tile to ncg columns
    
    # Calculate Cjj
    atm_wcss = []
    for cg in range(ncg):
        atmlabel = np.where( atm_cglabel == cg)[0]
        atm_wjj = 0.0
        for j1 in atmlabel:
            for j2 in atmlabel:
                atm_wjj += C[j1, j2] + lamb * D[j1, j2]
        atm_wjj /= atmlabel.size * atmlabel.size
        atm_wcss.append(atm_wjj)
    atm_wcss = np.tile(np.array(atm_wcss), (natm,1))
    #print(atm_wcss[0,:])
    
    # Calculate Cij and add to the atm_wcss matrix
    for cg in range(ncg):
        atmlabel = np.where( atm_cglabel == cg)[0]
        atm_wij = 0.0
        for jj in atmlabel:
            atm_wij += C[:, jj] + lamb * D[:, jj]
        atm_wij = atm_wij * 2 / atmlabel.size
        atm_wcss[:,cg] -= atm_wij
    #print(atm_wcss[0,:])
    
    # Calculate Cii and add to the atm_wcss matrix
    atm_wii = C.diagonal() + lamb * D.diagonal()
    atm_wii = np.tile(atm_wii.T, (ncg,1)).T
    atm_wcss += atm_wii
    #print(atm_wcss[0,:])
    
    return atm_wcss
    
def initialization_normalize(A, B, natm, ncg):
    # Iniitialization
    initial_list = []
    range_a = np.amax(A) - np.amin(A)
    range_b = np.amax(B) - np.amin(B)
    A = np.subtract(A, np.amin(A))
    B = np.subtract(B, np.amin(B))
    A = np.divide(A, range_a)
    B = np.divide(B, range_b)
    
    F = A + B
    for atm in range(natm):
        cor = F[atm,:]
        ind = [atm]
        for cg in range(ncg-1):
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
    
    return F, unique_initial_list
    
def convert_cglabel_to_atmlabel(ncg, atm_cglabel):
    cg_atmlabel = []
    for cg in range(ncg):
        atmlabel = np.where( atm_cglabel == cg)[0]
        atmlabel = np.add(atmlabel, 1)
        cg_atmlabel.append(atmlabel)
    
    return cg_atmlabel
    
def write_wcss_per_initial(filename, unique_initial_list, initial_wcss, ncg):
    unique_initial_list = np.add(unique_initial_list, 1)
    output_data = np.column_stack((unique_initial_list, np.array(initial_wcss)))
    fmt = ["%4d" for ii in range(ncg)]
    fmt.extend(["%10.6f"])
    np.savetxt(filename, output_data, fmt)
    
def write_site(filename, cg_atmlabel, ncg):
    print("kmeans> Here is the final CG sites:")
    with open(filename, "wb") as f:
        for nn, row in enumerate(cg_atmlabel):
            print(row)
            string = "Site_"+str(nn)+"="
            f.write(string.encode('ascii'))
            np.savetxt(f, [row], fmt='%4d')

# MAIN

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("ced", help="Input coved.dat")
    parser.add_argument("cod", help="Input covd.dat")
    parser.add_argument("natm", type=int, help="Number of atoms")
    parser.add_argument("ncg", type=int, help="Desired number of CG sites")

    args = parser.parse_args()
    # Decalre some frequently used variables
    natm = args.natm
    ncg = args.ncg
    print("kmeans> Number of atom (N): ", str(natm), " | Number of CG sites (K): ", str(ncg))

    # Read coved.dat
    print("kmeans> Reading " + str(args.ced))
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
    print("kmeans> Reading " + str(args.cod))
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


    # Iniitialization
    print("kmeans> Initialization starts")
    F, unique_initial_list = initialization_normalize(A, B, natm, ncg)
    print("kmeans> Initialization finished")
    print("kmeans> ", str(unique_initial_list.shape[0]), " initial sets found!")

    initial_wcss_bf = []
    initial_wcss_a_bf = []
    initial_wcss_b_bf = []
    # Sort the initial list
    for index_list in unique_initial_list:
        atm_wcss = F[:,index_list]
        atm_cglabel = np.argmin(atm_wcss, axis=1)

        # Before convergence, calculate residual WCSS
        cg_wcss_a = calc_residual_wcss(C, ncg, atm_cglabel)
        cg_wcss_b = calc_residual_wcss(D, ncg, atm_cglabel)
        cg_wcss_a_sum = cg_wcss_a.sum()
        cg_wcss_b_sum = cg_wcss_b.sum()
        lamb = cg_wcss_a_sum / cg_wcss_b_sum
        initial_wcss_a_bf.append(cg_wcss_a)
        initial_wcss_b_bf.append(cg_wcss_b)
        cg_wcss_b = np.multiply(cg_wcss_b, lamb)
        cg_wcss = np.add(cg_wcss_a, cg_wcss_b)
        initial_wcss_bf.append(cg_wcss.sum())
    index = [i[0] for i in sorted(enumerate(initial_wcss_bf), key=lambda x:x[1])]
    value = [i[1] for i in sorted(enumerate(initial_wcss_bf), key=lambda x:x[1])]
    print("kmeans> Here are top ten initials with the lowest WCSS")
    print(index[:10])
    sorted_unique_initial_list = unique_initial_list[index]

    print("kmeans> Write a full list of initial guess and final WCSS values")
    filename = "FWCSS_bf_cg" + str(ncg) + ".dat"
    write_wcss_per_initial(filename, sorted_unique_initial_list, sorted(initial_wcss_bf), ncg)


    # Update CG labels using WCSS
    # WCSS(C) = Cii - 2/S * Cij + 1/S^2 * Cjj
    # initial_wcss: { N_initial X 1 }
    initial_wcss_af = []
    initial_wcss_a_af = []
    initial_wcss_b_af = []
    initial_cglabel = []
    for kk, index_list in enumerate(sorted_unique_initial_list):
        lamb_iteration = []
        atm_wcss = F[:,index_list]
        
        # Declaring first atm_cglabel using initial guess
        atm_cglabel_new = np.argmin(atm_wcss, axis=1)
        atm_cglabel = np.zeros((natm,))
        if np.array_equal(atm_cglabel_new, atm_cglabel):
            print("kmeans> ERROR: Initial guess went wrong!")
            exit()
        else:
            atm_cglabel = atm_cglabel_new
            wcss_iteration = []
            wcss_iteration_a = []
            wcss_iteration_b = []
            cg_wcss_a = calc_residual_wcss(C, ncg, atm_cglabel)
            cg_wcss_b = calc_residual_wcss(D, ncg, atm_cglabel)
            cg_wcss_a_sum = cg_wcss_a.sum()
            cg_wcss_b_sum = cg_wcss_b.sum()
            lamb = cg_wcss_a_sum / cg_wcss_b_sum
            lamb_iteration.append(lamb)
            cg_wcss_b = np.multiply(cg_wcss_b, lamb)
            cg_wcss = np.add(cg_wcss_a, cg_wcss_b)
            wcss_iteration_a.append(cg_wcss_a_sum)
            wcss_iteration_b.append(cg_wcss_b_sum)
            wcss_iteration.append(cg_wcss.sum())
        
        # Main iteration: clustering according to WCSS per atom
        for jj in range(999):
            for ii in range(999):
                # Update atom's CG index
                # Break if converged
                atm_wcss = calc_wcss_per_atom(C, D, lamb, ncg, atm_cglabel)
                atm_cglabel_new = np.argmin(atm_wcss, axis=1)

                cg_wcss_a = calc_residual_wcss(C, ncg, atm_cglabel_new)
                cg_wcss_b = calc_residual_wcss(D, ncg, atm_cglabel_new)
                cg_wcss_a_sum = cg_wcss_a.sum()
                cg_wcss_b_sum = cg_wcss_b.sum()
                cg_wcss_b = np.multiply(cg_wcss_b, lamb)
                cg_wcss = np.add(cg_wcss_a, cg_wcss_b)
                wcss_iteration_a.append(cg_wcss_a_sum)
                wcss_iteration_b.append(cg_wcss_b_sum)
                wcss_iteration.append(cg_wcss.sum())

                if np.array_equal(atm_cglabel_new, atm_cglabel):
                    #print("Converged in ", str(ii), " iterations!")
                    # Save times of convergence
                    initial_converge_count.append(ii)
                    break
                else:
                    atm_cglabel = atm_cglabel_new
                    initial_converge_count = []
    
            filename = "FWCSS_cg" + str(ncg) + "_init" + str(kk+1) + "_iter" + str(jj+1) + ".dat"
            np.savetxt(filename, wcss_iteration, fmt="%1.6e")

            cg_wcss_a = calc_residual_wcss(C, ncg, atm_cglabel)
            cg_wcss_b = calc_residual_wcss(D, ncg, atm_cglabel)
            cg_wcss_a_sum = cg_wcss_a.sum()
            cg_wcss_b_sum = cg_wcss_b.sum()
            lamb_new = cg_wcss_a_sum / cg_wcss_b_sum
            if abs(lamb_new - lamb) / lamb * 100.0 < 5.0:
                filename = "AWCSS_cg" + str(ncg) + "_init" + str(kk+1) + "_iter.dat"
                np.savetxt(filename, wcss_iteration_a, fmt="%1.6e")
                filename = "BWCSS_cg" + str(ncg) + "_init" + str(kk+1) + "_iter.dat"
                np.savetxt(filename, wcss_iteration_b, fmt="%1.6e")
                filename = "FWCSS_cg" + str(ncg) + "_init" + str(kk+1) + "_iter.dat"
                np.savetxt(filename, wcss_iteration, fmt="%1.6e")
                break
            else:
                lamb = lamb_new
                lamb_iteration.append(lamb)
                wcss_iteration_a = []
                wcss_iteration_b = []
                wcss_iteration = []
        
        # After convergence, calculate residual WCSS (WCSS per CG set)
        initial_cglabel.append(atm_cglabel)
        cg_wcss_a = calc_residual_wcss(C, ncg, atm_cglabel)
        cg_wcss_b = calc_residual_wcss(D, ncg, atm_cglabel)
        initial_wcss_a_af.append(cg_wcss_a)
        initial_wcss_b_af.append(cg_wcss_b)
        cg_wcss_b = np.multiply(cg_wcss_b, lamb)
        cg_wcss = np.add(cg_wcss_a, cg_wcss_b)
        initial_wcss_af.append(cg_wcss.sum())
        filename = "lamb_cg" + str(ncg) + "_init" + str(kk+1) + "_iter.dat"
        np.savetxt(filename, lamb_iteration, fmt="%1.6e")

    print("kmeans> Maximum number of iteration is ", str(max(initial_converge_count)))
    print("kmeans> Write iteration counts")
    # Write iteration counts
    filename = "IterCount_cg" + str(ncg) + ".dat"
    np.savetxt(filename, initial_converge_count, fmt='%3d')
    print("kmeans> The smallest sum of residual WCSS is MINWCSS=", str(min(initial_wcss_af)))
    
    # Write WCSS_A per CG_site per initial_set
    print("kmeans> Write WCSS_A profile")
    filename = "AWCSS_bf_cg" + str(ncg) + ".dat"
    np.savetxt(filename, initial_wcss_a_bf, fmt="%1.6e")
    filename = "AWCSS_af_cg" + str(ncg) + ".dat"
    np.savetxt(filename, initial_wcss_a_af, fmt="%1.6e")
    
    # Write WCSS_B per CG_site per initial_set
    print("kmeans> Write WCSS_B profile")
    filename = "BWCSS_bf_cg" + str(ncg) + ".dat"
    np.savetxt(filename, initial_wcss_b_bf, fmt="%1.6e")
    filename = "BWCSS_af_cg" + str(ncg) + ".dat"
    np.savetxt(filename, initial_wcss_b_af, fmt="%1.6e")

    # Write WCSS per initial set
    print("kmeans> Write a full list of initial guess and final WCSS values")
    filename = "FWCSS_af_cg" + str(ncg) + ".dat"
    write_wcss_per_initial(filename, sorted_unique_initial_list, initial_wcss_af, ncg)
    
    # Find the best CG set (index) according to the smallest WCSS
    index = initial_wcss_af.index(min(initial_wcss_af))
    tm_cglabel = initial_cglabel[index]
    #print(atm_cglabel)

    # Convert cglabel to atmlabel
    cg_atmlabel = convert_cglabel_to_atmlabel(ncg, atm_cglabel)

    # Write sites
    filename = "SITE_cg" + str(ncg) + ".dat"
    print("kmeans> Write final CG sites to ", filename)
    write_site(filename, cg_atmlabel, ncg)

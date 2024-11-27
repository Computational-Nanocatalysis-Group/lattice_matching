#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 12:26:46 2021

@programming: Paulo C. D. Mendes
@revisions: Yizhen Song
@conceptualization and supervision: Sergey Kozlov

Please consider citing our article:
Yizhen Song, Paulo C. D. Mendes, Sergey Kozlov. J. Mater. Chem. A, 2023, 11, 13665-13676 
(https://doi.org/10.1039/d3ta01940c)


The function get_r_list() is inspired by ideas in: A. Zur and T. C. McGill, J. Appl. Phys., 1984, 55, 378–386
(https://doi.org/10.1063/1.333084)

"""

import numpy as np
import os
import sys


def get_r_list(area_s, area_f, max_area, tol=0.05, verbose = True):
    """
    returns a list of rs and rf values that satisfies:
    rs/rf = area_f/area_s with the constraints:
    rs <= Area_max/area_s and rf <= Area_max/area_f
    
    This function was written inspired by ideas in: A. Zur and T. C. McGill, J. Appl. Phys., 1984, 55, 378–386
    (https://doi.org/10.1063/1.333084). Please cite their study if you use this function.
    
    """    

    r_list = []
    rmax1 = int(max_area / area_s)
    rmax2 = int(max_area / area_f)
    if (verbose == True):
        print('rmax1, rmax2: ', rmax1, rmax2)
        print('area_f/area_s = ', area_f / area_s)

    approximation = 1000
    for rs in range(1, rmax1 + 1):
        for rf in range(1, rmax2 + 1):
            ratio = (float(rs)*area_s) / (float(rf)*area_f)
            if (ratio > 1-tol and ratio < 1+tol):
                r_list.append([rs, rf])
                if abs(float(rs) / float(rf) - area_f / area_s) < approximation:
                    rs_best,rf_best = rs, rf
                    approximation = abs(float(rs) / float(rf) - area_f / area_s)
    if not r_list:
        sys.exit('No acceptable rs/rf found within requested tolerance. Please consider increasing tol') 
    if (verbose == True):
        print("\nThe [rs, rf] pair that best approximates a2/a1 is: \n%d %d" %(rs_best,rf_best))
        print('Accuracy the best approximation: %f = rs/rf - a2/a1' %(approximation))
    return r_list, rs_best, rf_best


def get_areas(substrate_matrix, film_matrix):

    area_s = np.linalg.norm(np.cross(substrate_matrix[0],substrate_matrix[1]))
    area_f = np.linalg.norm(np.cross(film_matrix[0],film_matrix[1]))
    
    return area_s, area_f


def relative_deformations_method(atoms_substrate, atoms_film):
    supercell_s = atoms_substrate.copy()
    supercell_f = atoms_film.copy()
    
    Da = np.subtract(supercell_s[0], supercell_f[0])
    Db = np.subtract(supercell_s[1], supercell_f[1])
    Da_plus_Db = np.add(Da, Db)

    len_Da = np.linalg.norm(Da)/np.linalg.norm(supercell_s[0])
    len_Db = np.linalg.norm(Db)/np.linalg.norm(supercell_s[1])
    len_Da_plus_Db = np.linalg.norm(Da_plus_Db)/np.linalg.norm(np.add(supercell_s[0],supercell_s[1]))
    
    return len_Da, len_Db, len_Da_plus_Db


def scan_composite_combinations(substrate_matrix, film_matrix, Nsearch, tol, deformation_limit):
    area_s, area_f = get_areas(substrate_matrix, film_matrix)
    
    r_list, rs_best, rf_best = get_r_list(area_s,area_f, 
                        max_area=max(area_s*Nsearch**2,area_f*Nsearch**2), 
                        tol=tol, verbose = True)
    
    rs_list = []
    rf_list = []
    for rs, rf in r_list:
        rs_list.append(rs)
        rf_list.append(rf)

    T_substrate_rs_dict = {key: [] for key in rs_list}
    for a1s in range(-Nsearch,Nsearch+1):
        for a2s in range(-Nsearch,Nsearch+1):
            for b1s in range(-Nsearch,Nsearch+1):
                for b2s in range(-Nsearch,Nsearch+1):
                    if (a1s != 0):
                        det_sub = np.linalg.det([[a1s, a2s], [b1s, b2s]])                        
                        if (any(round(det_sub) == sublist[0] for sublist in r_list)):                            
                            T_substrate_rs_dict[round(det_sub)].append([a1s, a2s, b1s, b2s])
    T_film_rf_dict = {key: [] for key in rf_list}
    for a1f in range(-Nsearch,Nsearch+1):
        for a2f in range(-Nsearch,Nsearch+1):
            for b1f in range(-Nsearch,Nsearch+1):
                for b2f in range(-Nsearch,Nsearch+1):
                    if (a1f != 0):
                        det_film = np.linalg.det([[a1f, a2f], [b1f, b2f]])
                        if (any(round(det_film) == sublist[1] for sublist in r_list)):
                            T_film_rf_dict[round(det_film)].append([a1f, a2f, b1f, b2f])
                        
    min_deformation = 9999
    
    if os.path.exists("results_explored_matrices.txt"):
        os.remove("results_explored_matrices.txt")
    f = open("results_explored_matrices.txt", "a")
    f.write("a1s a2s b1s b2s a1f a2f b1f b2f |a|_s |b|_s |a|_f |b|_f area angle_ab max(Da,Db,Da+Db)" + "\n")

    for rsmatch,rfmatch in r_list:
        for a1s,a2s,b1s,b2s in T_substrate_rs_dict[rsmatch]:
            for a1f,a2f,b1f,b2f in T_film_rf_dict[rfmatch]:
                T_substrate = np.array([[a1s, a2s], [b1s, b2s]])
                substrate_supercell = np.matmul(T_substrate,substrate_matrix)
            
                T_film = np.array([[a1f, a2f], [b1f, b2f]])
                film_supercell = np.matmul(T_film, film_matrix)
                
                len_Da, len_Db, len_Da_plus_Db = relative_deformations_method(substrate_supercell, film_supercell)

                deformation = max(len_Da,len_Db,len_Da_plus_Db)
                                                                         
                if (deformation < min_deformation or deformation == min_deformation):
                    item_min_deformation = [substrate_supercell, film_supercell, deformation]
                    min_deformation_matrix = [a1s, a2s, b1s, b2s, a1f, a2f, b1f, b2f]
                    min_deformation = deformation
                if (deformation < deformation_limit):
                    len_a_s = np.linalg.norm(substrate_supercell[0])
                    len_b_s = np.linalg.norm(substrate_supercell[1])
                    len_a_f = np.linalg.norm(film_supercell[0])
                    len_b_f = np.linalg.norm(film_supercell[1])
                    area_final = np.linalg.norm(np.cross(substrate_supercell[0],substrate_supercell[1]))
                    angle_final = np.degrees(np.arccos(np.dot((substrate_supercell[0] / np.linalg.norm(substrate_supercell[0])), (substrate_supercell[1] / np.linalg.norm(substrate_supercell[1])))))
                    f.write("%d %d %d %d %d %d %d %d %.2f %.2f %.2f %.2f %d %d %f" %(a1s,a2s,b1s,b2s,a1f,a2f,b1f,b2f,len_a_s,len_b_s,len_a_f,len_b_f,area_final,angle_final,deformation) + "\n") 
    f.close()

    print()
    print("deformation = max(Da, Db, Da+Db)\n")
    print("min(deformation) = %f" %(item_min_deformation[2]))
    print("Was found applying:") 
    print("\nelements of the substrate transformation matrix: %s %s %s %s" %(min_deformation_matrix[0],min_deformation_matrix[1],min_deformation_matrix[2],min_deformation_matrix[3]))
    print("supercell lattice: \n%s" %(item_min_deformation[0]))
    print("elements of the film transformation matrix: %s %s %s %s" %(min_deformation_matrix[4],min_deformation_matrix[5],min_deformation_matrix[6],min_deformation_matrix[7]))
    print("  supercell lattice: \n%s" %(item_min_deformation[1]))
    print("\n\nFind all combinations that give this deformation in the file results_explored_matrices.txt")    
    print()

    return None


def main():

    #example of function call giving ZnO(0001) and Cu(111) matrices as surface unit cells where the match 4x3:3x3 is found.
    scan_composite_combinations(substrate_matrix=[[3.250, 0.000], [-1.625, 2.815]], film_matrix=[[2.553, 0.000], [1.276, 2.211]], Nsearch=4, tol=0.1,deformation_limit=0.1)
    
if __name__ == '__main__':
    main()
    

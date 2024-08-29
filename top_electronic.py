'''
This code extracts the atoms in the surface microstruture (regions 1, 2 and 3) of an adsorbate.
It will input the data for each of the atom and create a electronic_properties.csv file. The
electronic properties include the valence electrons, electronegativity and charge transfer of
each neighboring atoms.
'''

import os, shutil
from ase.visualize import view
from ase.io import read,write,xsf
from ase.geometry import get_distances
from ase.neighborlist import NeighborList,natural_cutoffs,build_neighbor_list,first_neighbors,get_connectivity_matrix,neighbor_list
import numpy as np
import pandas as pd
import csv

#define a function to create a list of folders of surfaces
def list_surface(directory):
    '''
    :param directory: main directory that contains the files
    :return folders: a list containg all the folders in the specified directory
    '''
    folders = [folder for folder in os.listdir(directory) if os.path.isdir(os.path.join(directory,folder))]
    try:
        folders.remove('pdos')
    except:
        print('no pdos folder')
    return folders

#define a function to create a list of folders of adsorption sites
def list_adsorbate(directory):
    '''
    :param directory: main directory that contains the files
    :return folders: a list containg all the folders in the specified directory
    '''
    folders = [folder for folder in os.listdir(directory) if os.path.isdir(os.path.join(directory,folder))]
    try:
        folders.remove('scf')
        folders.remove('dos')
        folders.remove('pp')
    except:
        print('no scf or pp')
    return folders

#define a function to compute the adorption energy
def eads(surface,adsorbate):
    '''
    param surface: surface name
    param adsorbate: adsorbate name
    return e_ads: adsorption energy
    '''
    e_co = -44.1124681422
    with open('{}/LOG'.format(surface),'r') as output_file:  # open OUTCAR located in the directories
        lines = reversed(output_file.readlines()) # read and store all lines into list from bottom to top line
        # search for final energy
        for line in lines:
            if 'Final energy' in line:
                line = line.split()[3]
                e_surface = float(line) # convert the string toten to float
                break
    with open('{}/{}/LOG'.format(surface,adsorbate),'r') as output_file:  # open OUTCAR located in the directories
        lines = reversed(output_file.readlines()) # read and store all lines into list from bottom to top line
        # search for final energy
        for line in lines:
            if 'Final energy' in line:
                line = line.split()[3]
                e_co_surface = float(line) # convert the string toten to float
                break
    e_ads = (e_co_surface - e_surface - e_co)*13.60
    return e_ads

#define a function to determine the atoms in local microstructures and their indices in the INPUT file
def neighbor(surface, adsorbate, min, max, atom_min, atom_max): #1.5 to 2.0 for r1
    '''
    param surface: surface name
    param adsorbate: adsorbate name
    param min: minimum distance
    param max: maximum distance
    return reg1: list of nearest neighboring atoms
    return reg1_i: list of indices of nearest neighboring atoms
    '''
    structure = read('{}/{}/INPUT-ads'.format(surface,adsorbate), format='espresso-in')
    atom_count=structure.get_chemical_symbols()
    for i in range(len(atom_count)):
        if(atom_count[i]=='C'):
            #print(i)
            index=i
            break
    reg = []
    reg_i = []
    for j in range(atom_min, atom_max):
        a = structure.get_distance(index,j,mic=True)
        if a > min and a < max:
            reg.append(structure.get_chemical_symbols()[j])
            reg_i.append(j+1)
    return reg,reg_i

#define a function to assign the valence electrons, electronegativity and bader charge of an atom in the monometallic surface
def eprop(r):
    '''
    param surface: name of surface
    param r: atom in a specified region
    return val: number of valence electrons of atom r
    return eneg: electronegativity of atom r
    return q_pure: bader charge of atom r in pure surface
    '''
    if r == "Co":
        val = 9
        eneg = 1.88
        q_pure = 17.01934167
    elif r == "Cu":
        val = 11
        eneg = 1.9
        q_pure = 19.015861
    elif r == "Fe":
        val = 8
        eneg  = 1.83
        q_pure = 16.01999511
    elif r == "Mn":
        val = 7
        eneg = 1.55
        q_pure = 15.03298589
    else:
        val = 10
        eneg = 1.91
        q_pure = 18.02582578
    return val,eneg,q_pure

# define a function to get the bader charge
def bader(surface, r_i, q_pure):
    '''
    param surface: name of surface
    param r_i: index of the aton
    param q_pure: bader charge of pure
    return q_net: bader charge density difference
    '''
    df_bader = pd.read_csv('../../bader-dos-anal-rev copy/{}/pp/ACF.dat'.format(surface), header=None, sep="\s+")
    df_bader = df_bader.drop(df_bader.index[[0,1]])
    df_bader = df_bader.reset_index(drop = True)
    q = df_bader.iloc[r_i-1,4]
    print(q)
    q_net = float(q) - q_pure
    print(q_net)
    return q_net


'''
start of main code
'''
# set main_path (directory where this py file is found). This py file is located in the same folder where the surface folders are located.
main_path = './'

# define a list of column names in a list.
data_hcp = [['surface','adsorbate','e_ads','r1_1','r1_2','r1_3','r2_1','r2_2','r2_3','r3_1','r1_1_i','r1_2_i','r1_3_i','r2_1_i','r2_2_i','r2_3_i','r3_1_i','val1_1','val1_2','val1_3','val2_1','val2_2','val2_3','val3_1','eneg1_1','eneg1_2','eneg1_3','eneg2_1','eneg2_2','eneg2_3','eneg3_1','q1_1','q1_2','q1_3','q2_1','q2_2','q2_3','q3_1']]

# open fcc and hcp csv files for directories names (surfaces and adsorbates)
ls_surface=list_surface(main_path)
for i in range(0,len(ls_surface)):
    surface = ls_surface[i]
    ls_adsorbate=list_adsorbate('{}/{}'.format(main_path,surface))
    for j in range(0,len(ls_adsorbate)):
        adsorbate = ls_adsorbate[j]
        print(surface, adsorbate)
        if 'x' in adsorbate or 'unf' in adsorbate or 'bridge' in adsorbate:
            print(surface,adsorbate,'data not included.')
        elif 'hcp' in adsorbate:
            try:
                e_ads= eads(surface,adsorbate)
                try:
                    r1, r1_i = neighbor(surface, adsorbate, 1.7, 2.0, 27, 36)
                    val1 = []
                    eneg1 = []
                    q_pure1 = []
                    q_net1 = []
                    for k in range(0,len(r1)):
                        val1.append(eprop(r1[k])[0])
                        eneg1.append(eprop(r1[k])[1])
                        q_pure1.append(eprop(r1[k])[2])
                        q_net1.append(bader(surface, r1_i[k], q_pure1[k]))
                    r2, r2_i = neighbor(surface, adsorbate,  2.1, 3.92, 27, 36)
                    val2 = []
                    eneg2 = []
                    q_pure2 = []
                    q_net2 = []
                    for k in range(0,len(r2)):
                        val2.append(eprop(r2[k])[0])
                        eneg2.append(eprop(r2[k])[1])
                        q_pure2.append(eprop(r2[k])[2])
                        q_net2.append(bader(surface, r2_i[k], q_pure2[k]))
                    r3, r3_i = neighbor(surface, adsorbate, 3.5, 4.54, 18, 27)
                    val3 = []
                    eneg3 = []
                    q_pure3 = []
                    q_net3 = []
                    for k in range(0,len(r3)):
                        val3.append(eprop(r3[k])[0])
                        eneg3.append(eprop(r3[k])[1])
                        q_pure3.append(eprop(r3[k])[2])
                        q_net3.append(bader(surface, r3_i[k], q_pure3[k]))
                    row = [surface,adsorbate,e_ads,r1[0],r2[0],r2[1],r2[2],r2[3],r2[4],r2[5],r3[0],r1_i[0],r1_i[1],r1_i[2],r2_i[0],r2_i[1],r2_i[2],r3_i[0],r3_i[1],r3_i[2],val1[0],val1[1],val1[2],val2[0],val2[1],val2[2],val3[0],eneg1[0],eneg1[1],eneg1[2],eneg2[0],eneg2[1],eneg2[2],eneg3[0],q_net1[0],q_net1[1],q_net1[2],q_net2[0],q_net2[1],q_net2[2],q_net3[0]]
                except:
                    print(surface+'_'+adsorbate+'does not have input-log.')
                    continue
            except:
                print(surface+'_'+adsorbate+'does not have input-log.')
                continue
            try:
                data_hcp.append(row)
            except:
                print('no row')
                continue

csv_file = 'electronic_hcp.csv' #change file name

#write the data to a CSV file
with open(csv_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(data_hcp)

print(f"CSV file '{csv_file}' has been created successfully.")

# sample
# surface = '25599029'
# adsorbate = 'fcc1'
# r1, r1_i = neighbor(surface, adsorbate, 1.5, 2.0, 27, 36)
# val = []
# eneg = []
# q_pure = []
# q_net = []
# for k in range(0,len(r1)):
#     print(eprop(r1[k]))
#     val.append(eprop(r1[k])[0])
#     eneg.append(eprop(r1[k])[1])
#     q_pure.append(eprop(r1[k])[2])
#     q_net.append(bader(surface, r1_i[k], q_pure[k]))
# row = [val[0],val[1],val[2],eneg[0],eneg[1],eneg[2],q_net[0],q_net[1],q_net[2]]
# print(row)


'''
end of code
'''

'''
This code extracts and compiles the PDOS (s, px, py, pz, dx2, dzx, dzy, dx2y2, and dxy) of the nearest neighboring atoms (regions 1, 2, and 3)
of an adsorbate in a surface. Each PDOS file can be found in <band>_dos.csv file inside the directory top_pdos/<region>.
'''


import os, shutil
from ase.visualize import view
from ase.io import read,write,xsf
from ase.geometry import get_distances
from ase.neighborlist import NeighborList,natural_cutoffs,build_neighbor_list,first_neighbors,get_connectivity_matrix,neighbor_list
import numpy as np
import pandas as pd
import csv


#define a function to create a list of folders
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

#define a function to create a list of folders
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

#define a function to get the fermi energy
def fermi(surface):
    '''
    :param main_path: main directory that contains the files
    :param surface : surface name
    :return fermi : fermi energy of the surface
    '''
    # extract fermi energy of the alloy from OUTCAR
    with open('../../bader-dos-anal-rev copy/{}/scf/LOG'.format(surface),'r') as output_file:  # open OUTCAR located in the directories
        lines = reversed(output_file.readlines()) # read and store all lines into list from bottom to top line
        # search for final energy
        for line in lines:
            if 'Fermi energy' in line:
                line = line.split()[4]
                fermi_fl = float(line) # convert the string toten to float
                break
    # returning energy of the alloy
    return fermi_fl

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
            index=i
            break
    reg = []
    reg_i = []
    reg_x = []
    reg_y = []
    reg_a = []
    for j in range(atom_min, atom_max):
        a = structure.get_distance(index,j,mic=True)
        if a > min and a < max:
            reg.append(structure.get_chemical_symbols()[j])
            reg_i.append(j+1)
            position = structure[j].position
            reg_x.append(position[0])
            reg_y.append(position[1])
            reg_a.append(a)

    return reg, reg_i, reg_x, reg_y, reg_a

def neighbor_clean(surface, r, reg_x, reg_y, start_index, end_index): #1.5 to 2.0 for r1
    structure_co = read('{}/{}/INPUT-ads'.format(surface,adsorbate), format='espresso-in')
    print(reg_x, reg_y)
    # Load your structure file
    try:
        structure = read('../../bader-dos-anal-rev copy/{}/scf/INPUT'.format(surface), format='espresso-in')
    except:
        structure = read('../../bader-dos-anal-rev copy/{}/scf/INPUT-ads'.format(surface), format='espresso-in')

    tol = 0.4

    for i in range(start_index, end_index):
        atom = structure[i]
        position = structure[i].position
        type = structure.get_chemical_symbols()[i]
        # print(atom)
        if reg_x-tol < position[0] and position[0] < reg_x+tol and reg_y-tol < position[1] and position[1] < reg_y+tol:
            if type == r:
                index = i
                # print('close:', type, i+1, position[0], position[1])
    return index+1

#define a function to extract the s, p, d pdos
def pdos(path, surface, adsorbate, region, fermi, r, r_i):
    '''
    :param surface: surface name
    :param adsorbate: adsorbate name
    :param fermi : fermi energy of the surface
    :param r1: atom
    :param r1_i: index of the atom
    '''
    min = -14
    max = 8
    #s-band
    df_s1 = pd.read_fwf('../../bader-dos-anal-rev copy/{}/dos/dos.pdos_atm#{}({})_wfc#1(s)'.format(surface, r_i, r))
    try:
        df_s2 = pd.read_fwf('../../bader-dos-anal-rev copy/{}/dos/dos.pdos_atm#{}({})_wfc#4(s)'.format(surface, r_i, r))
    except: #for Mn
        df_s2 = pd.read_fwf('../../bader-dos-anal-rev copy/{}/dos/dos.pdos_atm#{}({})_wfc#3(s)'.format(surface, r_i, r))
    df_s1['E-Ef'] = df_s1['# E (eV)'] - fermi
    df_s2['E-Ef'] = df_s2['# E (eV)'] - fermi
    df_s1_filtered = df_s1[(df_s1['E-Ef'] >= min) & (df_s1['E-Ef'] <= max)].copy()
    df_s2_filtered = df_s2[(df_s2['E-Ef'] >= min) & (df_s2['E-Ef'] <= max)].copy()
    #p-band
    df_p1 = pd.read_fwf('../../bader-dos-anal-rev copy/{}/dos/dos.pdos_atm#{}({})_wfc#2(p)'.format(surface, r_i, r))
    df_p2 = pd.read_fwf('../../bader-dos-anal-rev copy/{}/dos/dos.pdos_atm#{}({})_wfc#5(p)'.format(surface, r_i, r))
    df_p1['E-Ef'] = df_p1['# E (eV)'] - fermi
    df_p2['E-Ef'] = df_p2['# E (eV)'] - fermi
    df_p1_filtered = df_p1[(df_p1['E-Ef'] >= min) & (df_p1['E-Ef'] <= max)].copy()
    df_p2_filtered = df_p2[(df_p2['E-Ef'] >= min) & (df_p2['E-Ef'] <= max)].copy()
    #d-band
    try:
        df_d = pd.read_fwf('../../bader-dos-anal-rev copy/{}/dos/dos.pdos_atm#{}({})_wfc#3(d)'.format(surface, r_i, r))
    except: #for Mn
        df_d = pd.read_fwf('../../bader-dos-anal-rev copy/{}/dos/dos.pdos_atm#{}({})_wfc#4(d)'.format(surface, r_i, r))
    df_d['E-Ef'] = df_d['# E (eV)'] - fermi
    df_d_filtered = df_d[(df_d['E-Ef'] >= min) & (df_d['E-Ef'] <= max)].copy()
    df = df_s1_filtered[['E-Ef']].copy() #copy E-Ef column

    try:
        try:
            #adding all s pdos (s1_up,s1_dw,s2_up,s2_dw)
            df['s_pdos(E)'] = df_s1_filtered.iloc[:,3]+df_s2_filtered.iloc[:,3]+df_s1_filtered.iloc[:,4]+df_s2_filtered.iloc[:,4]
            #adding all p pdos (p1_up,p1_dw,p2_up,p2_dw)
            df['px_pdos(E)'] = df_p1_filtered.iloc[:,3]+df_p2_filtered.iloc[:,3]+df_p1_filtered.iloc[:,4]+df_p2_filtered.iloc[:,4]
            df['py_pdos(E)'] = df_p1_filtered.iloc[:,5]+df_p2_filtered.iloc[:,5]+df_p1_filtered.iloc[:,6]+df_p2_filtered.iloc[:,6]
            df['pz_pdos(E)'] = df_p1_filtered.iloc[:,7]+df_p2_filtered.iloc[:,7]+df_p1_filtered.iloc[:,8]+df_p2_filtered.iloc[:,8]
            #adding all d pdos (d1_up,d1_dw,d2_up,d2_dw)
            df['dz2_ldos(E)'] = df_d_filtered.iloc[:,3]+df_d_filtered.iloc[:,4]
            df['dzx_ldos(E)'] = df_d_filtered.iloc[:,5]+df_d_filtered.iloc[:,6]
            df['dzy_ldos(E)'] = df_d_filtered.iloc[:,7]+df_d_filtered.iloc[:,8]
            df['dx2y2_ldos(E)'] = df_d_filtered.iloc[:,9]+df_d_filtered.iloc[:,10]
            df['dxy_ldos(E)'] = df_d_filtered.iloc[:,11]+df_d_filtered.iloc[:,12]
        except:
            print(surface, adsorbate, region, fermi, r, r_i,'no dos')
    except:
        print(surface, adsorbate, region, fermi, r, r_i,'no dos')
    try:
        # add zero values of pdos to reach max E-Ef
        while df.loc[df.index[-1], 'E-Ef'] < max-0.01:
            new_rows = pd.DataFrame({'E-Ef': [df.loc[df.index[-1], 'E-Ef']+0.01],'s_pdos(E)': [0],
                'px_pdos(E)': [0],'py_pdos(E)': [0],'pz_pdos(E)': [0],
                'dz2_ldos(E)': [0],'dzx_ldos(E)': [0],'dzy_ldos(E)': [0],
                'dx2y2_ldos(E)': [0],'dxy_ldos(E)': [0]
            })
            df =pd.concat([df,new_rows], ignore_index=True)
    except:
        print(df.index[-1],'is the max value.')
    df.to_csv('{}/{}_{}_{}_{}_pdos.csv'.format(path,surface,adsorbate,r_i,r), index=False)


'''
start of main code
'''
#set main_path (directory where this py file is located). This .py file is located in the same folder where the surfaces folders are located.
main_path = './'
data_i = [['surface','adsorbate','r1_1', 'r2_1', 'r2_2', 'r2_3','r2_4', 'r2_5', 'r2_6', 'r3_1', 'r3_2', 'r3_3','r1_1_i', 'r2_1_i', 'r2_2_i', 'r2_3_i','r2_4_i', 'r2_5_i', 'r2_6_i', 'r3_1_i', 'r3_2_i', 'r3_3_i']]

# open fcc and top csv files for directories names (surfaces and adsorbates)
ls_surface=list_surface(main_path)
for i in range(0,len(ls_surface)):
    surface = ls_surface[i]
    ls_adsorbate=list_adsorbate('{}/{}'.format(main_path,surface))
    for j in range(0,len(ls_adsorbate)):
        adsorbate = ls_adsorbate[j]
        if adsorbate.startswith('top'):
            # print(surface, adsorbate)
            fermi_val = fermi(surface)
            reg1, reg1_i, reg1_x, reg1_y, reg1_a = neighbor(surface, adsorbate, 1.7, 2.0, 27, 36)
            reg1_value_pairs = list(zip(reg1, reg1_a, reg1_i, reg1_x, reg1_y))
            reg1_sorted = sorted(reg1_value_pairs, key=lambda x: x[1])
            r1 = [label for label, value1, value2, value3, value4 in reg1_sorted[:1]]
            r1_a = [value1 for label, value1, value2, value3, value4 in reg1_sorted[:1]]
            r1_i = [value2 for label, value1, value2, value3, value4 in reg1_sorted[:1]]
            r1_x = [value3 for label, value1, value2, value3, value4 in reg1_sorted[:1]]
            r1_y = [value4 for label, value1, value2, value3, value4 in reg1_sorted[:1]]
            for k in range(0,len(r1)):
                try:
                    pdos('../top_pdos/r1', surface, adsorbate, r1, fermi_val, r1[k], r1_i[k])
                except:
                    try:
                        print(r1[k], r1_i[k], r1_x[k], r1_y[k])
                        r1_i[k] = neighbor_clean(surface, r1[k], r1_x[k], r1_y[k], 27, 36)
                        print(r1[k], r1_i[k], r1_x[k], r1_y[k])
                        pdos('../top_pdos/r1', surface, adsorbate, r1, fermi_val, r1[k], r1_i[k])
                    except:
                        print(surface, adsorbate, 'skipped.')

            reg2, reg2_i, reg2_x, reg2_y, reg2_a = neighbor(surface, adsorbate,  2.1, 3.92, 27, 36)
            reg2_value_pairs = list(zip(reg2, reg2_a, reg2_i, reg2_x, reg2_y))
            reg2_sorted = sorted(reg2_value_pairs, key=lambda x: x[1])
            r2 = [label for label, value1, value2, value3, value4 in reg2_sorted[:6]]
            r2_a = [value1 for label, value1, value2, value3, value4 in reg2_sorted[:6]]
            r2_i = [value2 for label, value1, value2, value3, value4 in reg2_sorted[:6]]
            r2_x = [value3 for label, value1, value2, value3, value4 in reg2_sorted[:6]]
            r2_y = [value4 for label, value1, value2, value3, value4 in reg2_sorted[:6]]
            #print(r2)
            #print(r2_i)
            for k in range(0,len(r2)):
                try:
                    pdos('../top_pdos/r2', surface, adsorbate, r2, fermi_val, r2[k], r2_i[k])
                except:
                    try:
                        print(r2[k], r2_i[k], r2_x[k], r2_y[k])
                        r2_i[k] = neighbor_clean(surface, r2[k], r2_x[k], r2_y[k], 27, 36)
                        print(r2[k], r2_i[k], r2_x[k], r2_y[k])
                        pdos('../top_pdos/r2', surface, adsorbate, r2, fermi_val, r2[k], r2_i[k])
                    except:
                        print(surface, adsorbate, 'skipped.')

            reg3, reg3_i, reg3_x, reg3_y, reg3_a = neighbor(surface, adsorbate, 3.5, 4.54, 18, 27)
            reg3_value_pairs = list(zip(reg3, reg3_a, reg3_i, reg3_x, reg3_y))
            reg3_sorted = sorted(reg3_value_pairs, key=lambda x: x[1])
            r3 = [label for label, value1, value2, value3, value4 in reg3_sorted[:3]]
            r3_a = [value1 for label, value1, value2, value3, value4 in reg3_sorted[:3]]
            r3_i = [value2 for label, value1, value2, value3, value4 in reg3_sorted[:3]]
            r3_x = [value3 for label, value1, value2, value3, value4 in reg3_sorted[:3]]
            r3_y = [value4 for label, value1, value2, value3, value4 in reg3_sorted[:3]]
            #print(r3)
            #print(r3_i)
            for k in range(0,len(r3)):
                try:
                    pdos('../top_pdos/r3', surface, adsorbate, r3, fermi_val, r3[k], r3_i[k])
                except:
                    try:
                        print(r3[k], r2_i[k], r3_x[k], r3_y[k])
                        r3_i[k] = neighbor_clean(surface, r3[k], r3_x[k],r3_y[k], 18, 27)
                        print(r3[k], r3_i[k], r3_x[k], r3_y[k])
                        pdos('../top_pdos/r3', surface, adsorbate, r3, fermi_val, r3[k], r3_i[k])
                    except:
                        print(surface, adsorbate, 'skipped.')

            results_i = [str(ls_surface[i]), adsorbate, *r1, *r2, *r3, *r1_i, *r2_i, *r3_i]
            data_i.append(results_i)

i_file = 'electronic_top_hea.csv'
# Write the data to a CSV file
with open(i_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(data_i)

print(f"CSV file '{i_file}' has been created successfully.")

'''
end of code
'''

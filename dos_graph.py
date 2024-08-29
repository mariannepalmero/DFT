#PDOS plot

import matplotlib.pyplot as plt
import csv
import pandas as pd

# Lists to store data
x_values = []
s_up = []
s_dw = []
p_up = []
p_dw = []
d_up = []
d_dw = []

# Read data from CSV file
with open('fcc_pdos_pure/r1/25599029_fcc9_35_Ni_pdos.csv', 'r') as file: # change filename
    reader = csv.reader(file)
    next(reader)  # Skip header row if exists
    for row in reader:
        x_values.append(float(row[0]))  # first column contains x-values
        s_up.append(float(row[1]))
        s_dw.append(float(row[2]))
        p_up.append(float(row[3]))
        p_dw.append(float(row[4]))
        d_up.append(float(row[5]))
        d_dw.append(float(row[6]))
df = pd.read_csv('fcc_pdos_pure/r1/25599029_fcc9_35_Ni_pdos.csv') # change filename
df['total_dos_up'] = df['s_ldosup(E)']+ df['p_ldosup(E)']+ df['d_ldosup(E)']
df['total_dos_dw'] = df['s_ldosdw(E)']+ df['p_ldosdw(E)']+ df['d_ldosdw(E)']

# Plot the data

plt.plot(x_values, s_up, linestyle='-', label='s$_{up}$')
plt.plot(x_values, s_dw, linestyle='-', label='s$_{dw}$')
plt.plot(x_values, p_up, linestyle='-', label='p$_{up}$')
plt.plot(x_values, p_dw, linestyle='-', label='p$_{dw}$')
plt.plot(x_values, d_up, linestyle='-', label='d$_{up}$')
plt.plot(x_values, d_dw, linestyle='-', label='d$_{dw}$')
plt.plot(x_values, df['total_dos_up'], linestyle='-', label='total$_{up}$')
plt.plot(x_values, df['total_dos_dw'], linestyle='-', label='total$_{dw}$')
plt.ylabel('DOS')
plt.xlabel('E-E$_{f}$ (eV)')
plt.grid(True)
plt.ylim(-5, 5) # change range
plt.xlim(-8, 6) # change range
plt.legend(loc='upper left')
# plt.legend('s$_{up}$','s$_{down}$','p$_{up}$','p$_{down}$','d$_{up}$','d$_{down}$','total$_{up}$','total$_{down}$')
plt.savefig('pure_r1_25599029_fcc9_35_Ni_pdos.png') # change filename
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read data from CSV file using Pandas
df = pd.read_csv('complete-top-element.csv') # change the file name

# For different HEA system, change the element names of the following:
# Filter the data where 'Co1' is equal to 3
filter_co1 = df[df['Co1'] == 1]

# Filter the data where 'Cu1' is equal to 3
filter_cu1 = df[df['Cu1'] == 1]

# Filter the data where 'Fe1' is equal to 3
filter_fe1 = df[df['Fe1'] == 1]

# Filter the data where 'Mn1' is equal to 3
filter_mn1 = df[df['Mn1'] == 1]

# Filter the data where 'Ni1' is equal to 3
filter_ni1 = df[df['Ni1'] == 1]

# Extract the 'E_ads' column from the filtered data
e_ads_co1 = filter_co1['E_ads']
e_ads_cu1 = filter_cu1['E_ads']
e_ads_fe1 = filter_fe1['E_ads']
e_ads_mn1 = filter_mn1['E_ads']
e_ads_ni1 = filter_ni1['E_ads']
e_ads_all1 = df['E_ads']

# Calculate the bin edges based on the desired interval
bin_interval = 0.015  # Specify the interval for the bins
min_value = e_ads_all1.min()
max_value = e_ads_all1.max()
bin_edges = np.arange(min_value, max_value + bin_interval, bin_interval)

# Plot the histogram
plt.figure(figsize=(8, 6))
plt.hist(e_ads_all1, bins=bin_edges, color='cyan', alpha=0.7, edgecolor='black')
plt.hist(e_ads_co1, bins=bin_edges, color='blue', alpha=0.7, edgecolor='black')
plt.hist(e_ads_cu1, bins=bin_edges, color='green', alpha=0.7, edgecolor='black')
plt.hist(e_ads_fe1, bins=bin_edges, color='yellow', alpha=0.7, edgecolor='black')
plt.hist(e_ads_mn1, bins=bin_edges, color='magenta', alpha=0.7, edgecolor='black')
plt.hist(e_ads_ni1, bins=bin_edges, color='silver', alpha=0.7, edgecolor='black')

# For different HEA system, change the element names.
plt.legend(['all', '3Co', '3Cu', '3Fe', '3Mn', '3Ni'])

# Calculate the average of E_ads
ave_e_ads_all1 = np.mean(e_ads_all1)
ave_e_ads_co1 = np.mean(e_ads_co1)
ave_e_ads_cu1 = np.mean(e_ads_cu1)
ave_e_ads_fe1 = np.mean(e_ads_fe1)
ave_e_ads_mn1 = np.mean(e_ads_mn1)
ave_e_ads_ni1 = np.mean(e_ads_ni1)

print("Average E_ads Co:", ave_e_ads_co1)
print("Average E_ads Cu:", ave_e_ads_cu1)
print("Average E_ads Fe:", ave_e_ads_fe1)
print("Average E_ads Mn:", ave_e_ads_mn1)
print("Average E_ads Ni:", ave_e_ads_ni1)

# Plot a horizontal line for the average
# plt.axvline(x=ave_e_ads_all1, color='cyan', linestyle='--')
plt.axvline(x=ave_e_ads_co1, color='blue', linestyle='--')
plt.axvline(x=ave_e_ads_cu1, color='green', linestyle='--')
plt.axvline(x=ave_e_ads_fe1, color='yellow', linestyle='--')
plt.axvline(x=ave_e_ads_mn1, color='magenta', linestyle='--')
plt.axvline(x=ave_e_ads_ni1, color='silver', linestyle='--')
plt.axvline(x=0, color='black')

plt.xlabel('$E_{ads}$ (eV)', fontsize=18)
plt.ylabel('Frequency', fontsize=18)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.savefig('hist-top-element.jpg', dpi=300) # change the file name

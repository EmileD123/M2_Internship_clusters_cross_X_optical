import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import pickle
import os
import importlib
import sys
import pandas as pd


import halo_and_cluster_MF as hcmf
import useful_functions as uf

importlib.reload(hcmf) ; importlib.reload(uf)

# Le but de ce script est de tester les relations d'échelles fournies dans les différents papiers associés aux relevés utilisés, notamment pour vérifier si elles permettent de retrouver less mêmes valeurs

#----------------------------------------------------------------------------------------------------------------
# Test pour Gao et al. 2019 : SCUSS/SDSS/UNwise

'''here = os.getcwd()
print('here :', here)'''

# Valeurs

path_SCUSS = 'SCUSS_DESI/SCUSS_Gao_2020.fit'
table_SCUSS = uf.fits_to_dataframe(path_SCUSS)

zSCUSS_list = []
for idx, row in table_SCUSS.iterrows():
    if row[f'zsp'] == 0:
        zSCUSS_list.append(row[f'zph'])
    else:
        zSCUSS_list.append(row[f'zsp'])
table_SCUSS['z'] = zSCUSS_list

log10_M500_SCUSS_data = np.log10(np.array(table_SCUSS['M500']))
RL_SCUSS_data = np.array(table_SCUSS['Richness'])
L1Mpc_SCUSS_data = np.array(table_SCUSS['L1Mpc'])

log10_M500_SCUSS_retrieved_from_L1Mpc = np.array(uf.L1Mpc_to_log10_M_Msun(L1Mpc_SCUSS_data,zSCUSS_list))
log10_M500_SCUSS_retrieved_from_RL = np.array(uf.RL_to_log10_M_Msun(RL_SCUSS_data))

ratio_from_L1Mpc = log10_M500_SCUSS_data/log10_M500_SCUSS_retrieved_from_L1Mpc[0]

ratio_from_RL = (log10_M500_SCUSS_data/log10_M500_SCUSS_retrieved_from_RL)[0]
ratio_from_RL_min = (log10_M500_SCUSS_data/log10_M500_SCUSS_retrieved_from_RL)[1]
ratio_from_RL_max = (log10_M500_SCUSS_data/log10_M500_SCUSS_retrieved_from_RL)[2]

indexes_SCUSS = np.arange(0,len(table_SCUSS))

print(len(indexes_SCUSS))
print(np.shape(ratio_from_RL))

# Find indexes where ratios are outside the 0.9-1.1 range
outliers_L1Mpc_10percent = indexes_SCUSS[(ratio_from_L1Mpc > 1.1) | (ratio_from_L1Mpc < 0.9)]
outliers_RL_10percent = indexes_SCUSS[(ratio_from_RL > 1.1) | (ratio_from_RL < 0.9)]

outliers_L1Mpc_extreme = indexes_SCUSS[(ratio_from_L1Mpc > 1.9) | (ratio_from_L1Mpc < 0.1)]
outliers_RL_extreme = indexes_SCUSS[(ratio_from_RL > 1.9) | (ratio_from_RL < 0.1)]

table_SCUSS_outliers_RL_10percent = table_SCUSS.iloc[outliers_L1Mpc_10percent]

# Tracés / print

print("\nOutlier indexes for L1Mpc ratio (ratio > 1.1 or < 0.9):")
print(outliers_L1Mpc_10percent)
print("\nOutlier indexes for RL ratio (ratio > 1.1 or < 0.9):")
print(outliers_RL_10percent)
print("\nOutlier indexes for L1Mpc ratio (ratio > 1.9 or < 0.1):")
print(outliers_L1Mpc_extreme)
print("\nOutlier indexes for RL ratio (ratio > 1.9 or < 0.1):")
print(outliers_RL_extreme)

print('percentage of outliers L1Mpc 10% :', len(outliers_L1Mpc_10percent)/len(indexes_SCUSS)*100, '%')
print('percentage of outliers RL 10% :', len(outliers_RL_10percent)/len(indexes_SCUSS)*100, '%')
print('percentage of outliers L1Mpc extreme :', len(outliers_L1Mpc_extreme)/len(indexes_SCUSS)*100, '%')
print('percentage of outliers RL extreme :', len(outliers_RL_extreme)/len(indexes_SCUSS)*100, '%')

# Display the DataFrame in the terminal using pandas' to_string for better formatting
pd.set_option('display.max_rows', None)  # Show all rows (be careful with very large DataFrames)
pd.set_option('display.max_columns', None)  # Show all columns
pd.set_option('display.width', None)  # No line wrapping
print(table_SCUSS_outliers_RL_10percent.to_string(index=True))


plt.figure(figsize=(12, 6))
plt.scatter(indexes_SCUSS,ratio_from_RL,label='Data / retrieved',alpha=0.5)
#plt.scatter(indexes_SCUSS,ratio_from_RL_min,label='(Data / retrieved) min',alpha=0.5)
#plt.scatter(indexes_SCUSS,ratio_from_RL_max,label='(Data / retrieved) max',alpha=0.5)
plt.axhline(y=1,color='red',label='1')
plt.axhline(y=0.9,color='orange',label='0.9')
plt.axhline(y=1.1,color='green',label='1.1')
plt.xlabel('Index')
plt.ylabel('RL data / RL retrieved')
plt.show()

plt.figure(figsize=(12, 6))
plt.scatter(indexes_SCUSS,ratio_from_L1Mpc,label='Data / retrieved',alpha=0.5)
plt.axhline(y=1,color='red',label='1')
plt.axhline(y=0.9,color='orange',label='0.9')
plt.axhline(y=1.1,color='green',label='1.1')
plt.xlabel('Index')
plt.ylabel('L1Mpc data / L1Mpc retrieved')
plt.show()
#----------------------------------------------------------------------------------------------------------------











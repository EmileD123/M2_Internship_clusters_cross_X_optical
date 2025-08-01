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

# Le but de ce script est de tester les relations d'échelles fournies dans les différents papiers associés aux relevés utilisés, notamment pour vérifier si elles permettent de retrouver les mêmes valeurs

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

log10_M500_SCUSS_data = np.log10(np.array(table_SCUSS['M500'])) ; M500_SCUSS_data = np.array(table_SCUSS['M500'])
RL_SCUSS_data = np.array(table_SCUSS['Richness'])
L1Mpc_SCUSS_data = np.array(table_SCUSS['L1Mpc'])

log10_M500_SCUSS_retrieved_from_L1Mpc = np.array(uf.L1Mpc_to_log10_M_Msun_test(L1Mpc_SCUSS_data,zSCUSS_list)) ; M500_SCUSS_retrieved_from_L1Mpc = 10**log10_M500_SCUSS_retrieved_from_L1Mpc
log10_M500_SCUSS_retrieved_from_RL = np.array(uf.RL_to_log10_M_Msun_test(RL_SCUSS_data)) ; M500_SCUSS_retrieved_from_RL = 10**log10_M500_SCUSS_retrieved_from_RL

ratio_log10_from_L1Mpc = log10_M500_SCUSS_data/log10_M500_SCUSS_retrieved_from_L1Mpc ; ratio_diff_from_L1Mpc = ((M500_SCUSS_data-M500_SCUSS_retrieved_from_L1Mpc)/M500_SCUSS_data)*100

ratio_log10_from_RL = (log10_M500_SCUSS_data/log10_M500_SCUSS_retrieved_from_RL)[0] ; ratio_diff_from_RL = ((M500_SCUSS_data-M500_SCUSS_retrieved_from_RL)/M500_SCUSS_data)[0]*100
ratio_log10_from_RL_min = (log10_M500_SCUSS_data/log10_M500_SCUSS_retrieved_from_RL)[1] ; ratio_from_RL_min = ((M500_SCUSS_data-M500_SCUSS_retrieved_from_RL)/M500_SCUSS_data)[1]*100
ratio_log10_from_RL_max = (log10_M500_SCUSS_data/log10_M500_SCUSS_retrieved_from_RL)[2] ; ratio_from_RL_max = ((M500_SCUSS_data-M500_SCUSS_retrieved_from_RL)/M500_SCUSS_data)[2]*100

indexes_SCUSS = np.arange(0,len(table_SCUSS))

# Find indexes where ratios are outside the 0.9-1.1 range
outliers_L1Mpc_10percent = indexes_SCUSS[(ratio_log10_from_L1Mpc > 1.1) | (ratio_log10_from_L1Mpc < 0.9)]
outliers_RL_10percent = indexes_SCUSS[(ratio_log10_from_RL > 1.1) | (ratio_log10_from_RL < 0.9)]

outliers_L1Mpc_extreme = indexes_SCUSS[(ratio_log10_from_L1Mpc > 1.9) | (ratio_log10_from_L1Mpc < 0.1)]
outliers_RL_extreme = indexes_SCUSS[(ratio_log10_from_RL > 1.9) | (ratio_log10_from_RL < 0.1)]

outliers_L1Mpc_null = indexes_SCUSS[ratio_log10_from_L1Mpc == 0]
outliers_RL_null = indexes_SCUSS[ratio_log10_from_RL == 0]


table_SCUSS_outliers_RL_10percent = table_SCUSS.iloc[outliers_L1Mpc_10percent]
table_SCUSS_outliers_RL_extreme = table_SCUSS.iloc[outliers_RL_extreme]
table_SCUSS_outliers_RL_null = table_SCUSS.iloc[outliers_RL_null]






# Tracés / print

print('Total size of ratio_log10_from_RL:', ratio_log10_from_RL.size)
print('Total size of ratio_log10_from_L1Mpc :', ratio_log10_from_L1Mpc.size)

print("Stats for ratio_from_L1Mpc:")
print(f"  Mean: {np.mean(ratio_diff_from_L1Mpc)}")
print(f"  Median: {np.median(ratio_diff_from_L1Mpc)}")
print(f"  Std: {np.std(ratio_diff_from_L1Mpc)}")
print(f"  Min: {np.min(ratio_diff_from_L1Mpc)}")
print(f"  Max: {np.max(ratio_diff_from_L1Mpc)}")
print(f"  25th percentile: {np.percentile(ratio_diff_from_L1Mpc, 25)}")
print(f"  75th percentile: {np.percentile(ratio_diff_from_L1Mpc, 75)}")

print("Stats for ratio_from_RL:")
print(f"  Mean: {np.mean(ratio_diff_from_RL)}")
print(f"  Median: {np.median(ratio_diff_from_RL)}")
print(f"  Std: {np.std(ratio_diff_from_RL)}")
print(f"  Min: {np.min(ratio_diff_from_RL)}")
print(f"  Max: {np.max(ratio_diff_from_RL)}")
print(f"  25th percentile: {np.percentile(ratio_diff_from_RL, 25)}")
print(f"  75th percentile: {np.percentile(ratio_diff_from_RL, 75)}")

'''print("\nOutlier indexes for L1Mpc ratio_log10 (ratio_log10 > 1.1 or < 0.9):")
print(outliers_L1Mpc_10percent)
print("\nOutlier indexes for RL ratio_log10 (ratio_log10 > 1.1 or < 0.9):")
print(outliers_RL_10percent)
print("\nOutlier indexes for L1Mpc ratio_log10 (ratio_log10 > 1.9 or < 0.1):")
print(outliers_L1Mpc_extreme)
print("\nOutlier indexes for RL ratio_log10 (ratio_log10 > 1.9 or < 0.1):")
print(outliers_RL_extreme)
print("\nOutlier indexes for L1Mpc ratio_log10 (ratio_log10 == 0):")
print(outliers_L1Mpc_null)
print("\nOutlier indexes for RL ratio_log10 (ratio_log10 == 0):")
print(outliers_RL_null)'''

'''print('percentage of outliers L1Mpc 10% :', len(outliers_L1Mpc_10percent)/len(indexes_SCUSS)*100, '%')
print('percentage of outliers RL 10% :', len(outliers_RL_10percent)/len(indexes_SCUSS)*100, '%')
print('percentage of outliers L1Mpc extreme :', len(outliers_L1Mpc_extreme)/len(indexes_SCUSS)*100, '%')
print('percentage of outliers RL extreme :', len(outliers_RL_extreme)/len(indexes_SCUSS)*100, '%')
print('percentage of outliers L1Mpc null :', len(outliers_L1Mpc_null)/len(indexes_SCUSS)*100, '%')
print('percentage of outliers RL null :', len(outliers_RL_null)/len(indexes_SCUSS)*100, '%')'''

# Display the DataFrame in the terminal using pandas' to_string for better formatting
pd.set_option('display.max_rows', None)  # Show all rows (be careful with very large DataFrames)
pd.set_option('display.max_columns', None)  # Show all columns
pd.set_option('display.width', None)  # No line wrapping
print("table_SCUSS_outliers_RL_null :", table_SCUSS_outliers_RL_null.to_string(index=True))

print('')


'''
# log10
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# First subplot: RL
axes[0].scatter(indexes_SCUSS, ratio_log10_from_RL, label='Data / retrieved', alpha=0.5)
#axes[0].scatter(indexes_SCUSS, ratio_log10_from_RL_min, label='(Data / retrieved) min', alpha=0.5)
#axes[0].scatter(indexes_SCUSS, ratio_log10_from_RL_max, label='(Data / retrieved) max', alpha=0.5)
axes[0].axhline(y=1, color='red', label='1')
axes[0].axhline(y=0.9, color='orange', label='0.9')
axes[0].axhline(y=1.1, color='green', label='1.1')
axes[0].set_xlabel('Index')
axes[0].set_ylabel('RL data / RL retrieved')
axes[0].legend()
axes[0].set_title('RL ratio_log10')

# Second subplot: L1Mpc
axes[1].scatter(indexes_SCUSS, ratio_log10_from_L1Mpc, label='Data / retrieved', alpha=0.5)
axes[1].axhline(y=1, color='red', label='1')
axes[1].axhline(y=0.9, color='orange', label='0.9')
axes[1].axhline(y=1.1, color='green', label='1.1')
axes[1].set_xlabel('Index')
axes[1].set_ylabel('L1Mpc data / L1Mpc retrieved')
axes[1].legend()
axes[1].set_title('L1Mpc ratio_log10')

plt.tight_layout()
plt.show()'''

# idem sans log10
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# First subplot: RL
axes[0].scatter(indexes_SCUSS, ratio_diff_from_RL, label='(Data - retrieved) / Data  (%)', alpha=0.5)
#axes[0].scatter(indexes_SCUSS, ratio_from_L1Mpc, label='Data / retrieved L1Mpc', alpha=0.5)
#axes[0].scatter(indexes_SCUSS, ratio_from_RL_min, label='(Data / retrieved) min', alpha=0.5)
#axes[0].scatter(indexes_SCUSS, ratio_from_RL_max, label='(Data / retrieved) max', alpha=0.5)
axes[0].axhline(y=0 , color='red')
axes[0].axhline(y = 5, color='orange')
axes[0].axhline(y=-5, color='green')
axes[0].set_xlabel('Index')
axes[0].set_ylabel('(RL data - RL retrieved) / RL data (%)')
axes[0].legend()
axes[0].set_title('RL ratio')

# Second subplot: L1Mpc
axes[1].scatter(indexes_SCUSS, ratio_diff_from_L1Mpc, label='(Data - retrieved) / Data (%)', alpha=0.5)
axes[1].axhline(y=0 , color='red')
axes[1].axhline(y = 5, color='orange')
axes[1].axhline(y=-5, color='green')
axes[1].set_xlabel('Index')
axes[1].set_ylabel('(L1Mpc data - L1Mpc retrieved) / L1Mpc data (%)')
axes[1].legend()
axes[1].set_title('L1Mpc ratio')

plt.tight_layout()
plt.show()
#----------------------------------------------------------------------------------------------------------------











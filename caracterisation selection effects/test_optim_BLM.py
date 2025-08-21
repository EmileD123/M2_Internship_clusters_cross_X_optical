import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import pickle
import os
import importlib
import sys
from astropy.cosmology import Planck18



import halo_and_cluster_MF as hcmf
import useful_functions as uf
from scipy.optimize import minimize

importlib.reload(hcmf) ; importlib.reload(uf)

cosmo = Planck18

path_eRASS1 = './DESI_eRass/eRass_Bulbul_2024_all.fit' #la première table est la table eRosita de départ, enrichie avec des informations d'association avec SDSS-Wen
table_eRASS1 = uf.fits_to_dataframe(path_eRASS1)

table_eRASS1 = table_eRASS1[table_eRASS1['F500'] < 10000]
table_eRASS1 = table_eRASS1[table_eRASS1['M500'] > 1]


M500_eRASS1_for_optim = np.array(table_eRASS1['M500']*1e13)
z_eRASS1_for_optim = np.array(table_eRASS1['zBest'])
Lx_02_23_for_optim = np.array(table_eRASS1['L500'])#*1e42)

plt.figure(figsize=(8, 6))
L_Maughan = (uf.L_from_M_Maughan_test_B_LM(M500_eRASS1_for_optim, z_eRASS1_for_optim, B_LM=1.64)*(1/2.57)) / 1e42
plt.hist(L_Maughan/Lx_02_23_for_optim - 1, bins=50, range=[-1,1])
plt.title('Histogram of L_Maughan / Lx_02_23 sans optim')
plt.legend()
#plt.show()





def cost_fun(params):
    B_LM, coeff = params
    L_Maughan = (uf.L_from_M_Maughan_test_B_LM(M500_eRASS1_for_optim, z_eRASS1_for_optim, B_LM=B_LM) * coeff) / 1e42
    return np.sqrt(np.sum((L_Maughan/Lx_02_23_for_optim - 1)**2))

def cost_fun_diff(params):
    B_LM, coeff = params
    L_Maughan = (uf.L_from_M_Maughan_test_B_LM(M500_eRASS1_for_optim, z_eRASS1_for_optim, B_LM=B_LM) * coeff) / 1e42
    return np.sum(((L_Maughan-Lx_02_23_for_optim)/(0.1*L_Maughan))**2)

result = minimize(cost_fun, x0=[1.3, 1/2])
B_LM_optim, coeff_optim = result.x

print(f'Optimized B_LM = {B_LM_optim}')
print(f'Optimized Coeff = {coeff_optim}')

L_Maughan = (uf.L_from_M_Maughan_test_B_LM(M500_eRASS1_for_optim, z_eRASS1_for_optim, B_LM=B_LM_optim)*(coeff_optim)) / 1e42

plt.figure(figsize=(8, 6))
plt.hist(L_Maughan/Lx_02_23_for_optim - 1, bins=50, range=[-1,1])
plt.title('Histogram of L_Maughan / Lx_02_23 pour optim 1')
plt.legend()
#plt.show()




ratio_L_Maughan_Lx_02_23 = L_Maughan/Lx_02_23_for_optim

table_eRASS1['ratio_L'] = ratio_L_Maughan_Lx_02_23

#plt.hist(ratio_L_Maughan_Lx_02_23,bins=50, range=[-0.1,2.1])

table_eRASS1_filtered = table_eRASS1[~np.isnan(table_eRASS1['ratio_L'])]
ratio_L_Maughan_Lx_02_23_filtered = table_eRASS1_filtered['ratio_L'] 


M500_eRASS1_for_optim = np.array(table_eRASS1_filtered['M500']*1e13)
z_eRASS1_for_optim = np.array(table_eRASS1_filtered['zBest'])
Lx_02_23_for_optim = np.array(table_eRASS1_filtered['L500'])#*1e42)


print("test cost_fun 2 = ",cost_fun([1.3, 1/2]))
result_2 = minimize(cost_fun, x0=[1.3, 1/3])
B_LM_optim_2, coeff_optim_2 = result_2.x
L_Maughan = (uf.L_from_M_Maughan_test_B_LM(M500_eRASS1_for_optim, z_eRASS1_for_optim, B_LM=B_LM_optim_2)*(coeff_optim_2)) / 1e42

print(f'---Optim result {result_2}')
print(f'Optimized B_LM 2 = {B_LM_optim_2}')
print(f'Optimized Coeff 2= {coeff_optim_2}')

M500_eRASS1_for_optim = np.array(table_eRASS1['M500']*1e13)
z_eRASS1_for_optim = np.array(table_eRASS1['zBest'])
Lx_02_23_for_optim = np.array(table_eRASS1['L500'])#*1e42)


L_Maughan = (uf.L_from_M_Maughan_test_B_LM(M500_eRASS1_for_optim, z_eRASS1_for_optim, B_LM=B_LM_optim_2)*(coeff_optim_2)) / 1e42
#plt.scatter(Lx_02_23_for_optim, diff_L_Maughan_Lx_02_23, alpha=0.1)
plt.figure(figsize=(8, 6))
plt.hist(L_Maughan/Lx_02_23_for_optim - 1, bins=50, range=[-1,1])
plt.title('Histogram of L_Maughan / Lx_02_23 pour optim 2')
plt.legend()
#plt.show()







diff_L_Maughan_Lx_02_23 = np.abs((L_Maughan-Lx_02_23_for_optim)/(0.1*L_Maughan))

table_eRASS1['diff_L'] = diff_L_Maughan_Lx_02_23
table_eRASS1_filtered_diff = table_eRASS1[table_eRASS1['diff_L'] < 50]

M500_eRASS1_for_optim = np.array(table_eRASS1_filtered_diff['M500']*1e13)
z_eRASS1_for_optim = np.array(table_eRASS1_filtered_diff['zBest'])
Lx_02_23_for_optim = np.array(table_eRASS1_filtered_diff['L500'])#*1e42)

print('sortie appel cost_fun_diff')
result_3 = minimize(cost_fun_diff, x0=[1.3, 1/3])
B_LM_optim_3, coeff_optim_3 = result_3.x
print(f'---Optim result diff {result_3}')
print(f'Optimized B_LM 3 = {B_LM_optim_3}')
print(f'Optimized Coeff 3= {coeff_optim_3}')

L_Maughan = (uf.L_from_M_Maughan_test_B_LM(M500_eRASS1_for_optim, z_eRASS1_for_optim, B_LM=B_LM_optim_3)*(coeff_optim_3)) / 1e42

#plt.scatter(Lx_02_23_for_optim,(L_Maughan-Lx_02_23_for_optim)/(0.1*L_Maughan),alpha=0.1)
plt.figure(figsize=(8, 6))
plt.hist(L_Maughan/Lx_02_23_for_optim - 1, bins=50, range=[-1,1])
plt.title('Histogram of L_Maughan / Lx_02_23 pour optim 3')
plt.legend()
plt.show()



L_Maughan_filtered = uf.L_from_M_Maughan_test_B_LM(table_eRASS1_filtered['M500']*1e13, table_eRASS1_filtered['zBest'], B_LM=B_LM_optim_2) * coeff_optim_2
ratio_L_Maughan_Lx_02_23_filtered_2 = L_Maughan_filtered/(table_eRASS1_filtered['L500'] * 1e42)

print(f'Mean ratio L_Maughan / Lx_02_23 = {np.mean(ratio_L_Maughan_Lx_02_23_filtered_2)}')
print(f'Standard deviation ratio L_Maughan / Lx_02_23 = {np.std(ratio_L_Maughan_Lx_02_23_filtered_2)}')



'''L_Maughan_filtered = uf.L_from_M_Maughan_test_B_LM(table_eRASS1_filtered['M500']*1e13, table_eRASS1_filtered['zBest'], B_LM=1.6) * 1/2.5
ratio_L_Maughan_Lx_02_23_filtered_2 = L_Maughan_filtered/(table_eRASS1_filtered['L500'] * 1e42)'''



'''plt.figure(figsize=(8, 6))
plt.scatter(table_eRASS1_filtered['L500'] * 1e42, ratio_L_Maughan_Lx_02_23_filtered_2, alpha=0.1, label='L_Maughan / Lx_02_23')
plt.axhline(y=1, color='red')
plt.axhline(y=1.05, color='orange')
plt.axhline(y=0.95, color='green')
plt.xlabel('Lx_02_23 [erg/s]')
plt.ylabel('L_Maughan / Lx_02_23 filtered')
plt.xscale('log')
plt.legend()
plt.show()

print(f'Mean ratio L_Maughan / Lx_02_23 = {np.mean(ratio_L_Maughan_Lx_02_23_filtered)}')
print(f'Standard deviation ratio L_Maughan / Lx_02_23 = {np.std(ratio_L_Maughan_Lx_02_23_filtered)}')

plt.figure(figsize=(8, 6))
plt.scatter(Lx_02_23, ratio_L_Maughan_Lx_02_23, c=z_eRASS1, cmap='viridis', alpha=0.1, label='L_Maughan / Lx_02_23')
plt.colorbar(label='Redshift z')
plt.axhline(y=1, color='red')
plt.axhline(y=1.05, color='orange')
plt.axhline(y=0.95, color='green')
plt.xlabel('Lx_02_23 [erg/s]')
plt.ylabel('L_Maughan / Lx_02_23')
plt.xscale('log')
plt.legend()
plt.show()'''

# Le but du script est d'optimiser la valeur du paramètre B_LM présent dans la scaling relation L-M du papier Maughan et al. 2018 ainsi que la valeur du coefficient permettant de passer à une luminosité dans 0.2-2.3 kev
# Pour ce faire, on va créer des fonctions coûts comparant les résultats de la scaling relation aux données en Lx de Bulbul et al. 2024

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

path_eRASS1 = './DESI_eRass/eRass_Bulbul_2024.fit' #la première table est la table eRosita de départ, enrichie avec des informations d'association avec SDSS-Wen
table_eRASS1 = uf.fits_to_dataframe(path_eRASS1)

# On commence par filtrer la table
table_eRASS1 = table_eRASS1[table_eRASS1['F500'] < 10000]
table_eRASS1 = table_eRASS1[table_eRASS1['M500'] > 1e-9] # 1805 elements of the table have a total mass M500 null ! → for them the selection funcion is "not reliable"
                                                         # More precisely that means that the scaling relations they used couldn't apply to these clusters (not enough counts etc ...)


# On récupère 2 sous-tables de eRASS1 avec des valeurs de z croissantes, contenant toutes les deux à peu près le même nombre d'élements
[table_eRASS1_lowz, table_eRASS1_highz] = uf.split_dataframe_in_k_bins_z(table_eRASS1, 'zBest', 2)
# ↑ Ca va nous permettre de définir une évolution linaire en fonction de z de B_LM et du z !

print(('Mean z in low_z = ', np.mean(table_eRASS1_lowz['zBest'])))
print(('Mean z in high_z = ', np.mean(table_eRASS1_highz['zBest'])))


# On récupère les données utiles
M500_eRASS1_for_optim = np.array(table_eRASS1['M500']*1e13)
z_eRASS1_for_optim = np.array(table_eRASS1['zBest'])
Lx_02_23_for_optim = np.array(table_eRASS1['L500'])#*1e42)

M500_eRASS1_for_optim_low_z = np.array(table_eRASS1_lowz['M500']*1e13)
z_eRASS1_for_optim_low_z = np.array(table_eRASS1_lowz['zBest'])
Lx_02_23_for_optim_low_z = np.array(table_eRASS1_lowz['L500'])#*1e42)

M500_eRASS1_for_optim_high_z = np.array(table_eRASS1_highz['M500']*1e13)
z_eRASS1_for_optim_high_z = np.array(table_eRASS1_highz['zBest'])
Lx_02_23_for_optim_high_z = np.array(table_eRASS1_highz['L500'])#*1e42)




# On calcule et on trace pour la luminosité avec le B_LM de base du papier et le coeff le plus adapté dans ce cas
plt.figure(figsize=(8, 6))
L_Maughan_sans_optim = (uf.L_from_M_Maughan_test_B_LM(M500_eRASS1_for_optim, z_eRASS1_for_optim, B_LM=1.64)*(1/2.57)) / 1e42

ratio_L_Maughan_Lx_02_23_sans_optim = L_Maughan_sans_optim/Lx_02_23_for_optim

print(f'Mean ratio L_Maughan / Lx_02_23 sans optim= {np.mean(ratio_L_Maughan_Lx_02_23_sans_optim)}')
print(f'Standard deviation ratio L_Maughan / Lx_02_23 sans optim= {np.std(ratio_L_Maughan_Lx_02_23_sans_optim)}')

plt.hist(ratio_L_Maughan_Lx_02_23_sans_optim - 1, bins=50, range=[-1,1])
plt.title('Histogram of L_Maughan / Lx_02_23 sans optim')
plt.legend()
#plt.show()




# On définit les fonctions coût
def cost_fun(params):
    B_LM, coeff = params
    L_Maughan = (uf.L_from_M_Maughan_test_B_LM(M500_eRASS1_for_optim, z_eRASS1_for_optim, B_LM=B_LM) * coeff) / 1e42 # On divise par 1e42 pour que ce soit comparable à Lx_02_23_for_optim : on reste dans cette ordre de grandeur pour faciliter l'optim
    return np.sqrt(np.sum((L_Maughan/Lx_02_23_for_optim - 1)**2))

def cost_fun_diff(params):
    B_LM, coeff = params
    L_Maughan = (uf.L_from_M_Maughan_test_B_LM(M500_eRASS1_for_optim, z_eRASS1_for_optim, B_LM=B_LM) * coeff) / 1e42
    return np.sum(((L_Maughan-Lx_02_23_for_optim)/(0.1*L_Maughan))**2)





# Première optimisation : cost_fun pour table entière
result = minimize(cost_fun, x0=[1.6, 1/2.5])
B_LM_optim, coeff_optim = result.x

print(f'---Optim 1 result {result}')
print(f'Optimized B_LM = {B_LM_optim}')
print(f'Optimized Coeff = {coeff_optim}')

L_Maughan = (uf.L_from_M_Maughan_test_B_LM(M500_eRASS1_for_optim, z_eRASS1_for_optim, B_LM=B_LM_optim)*(coeff_optim)) / 1e42
ratio_L_Maughan_Lx_02_23 = L_Maughan/Lx_02_23_for_optim
table_eRASS1['ratio_L'] = ratio_L_Maughan_Lx_02_23

print(f'Mean ratio L_Maughan / Lx_02_23 = {np.mean(ratio_L_Maughan_Lx_02_23)}')
print(f'Standard deviation ratio L_Maughan / Lx_02_23 = {np.std(ratio_L_Maughan_Lx_02_23)}')

plt.figure(figsize=(8, 6))
plt.hist(ratio_L_Maughan_Lx_02_23 - 1, bins=50, range=[-1,1])
plt.title('Histogram of L_Maughan / Lx_02_23 pour optim 1')
plt.legend()
#plt.show()





# Seconde optimisation : cost_fun_diff pour table entière
diff_L_Maughan_Lx_02_23 = np.abs((L_Maughan-Lx_02_23_for_optim)/(0.1*L_Maughan))

table_eRASS1['diff_L'] = diff_L_Maughan_Lx_02_23
table_eRASS1_filtered_diff = table_eRASS1[table_eRASS1['diff_L'] < 50]

M500_eRASS1_for_optim = np.array(table_eRASS1_filtered_diff['M500']*1e13)
z_eRASS1_for_optim = np.array(table_eRASS1_filtered_diff['zBest'])
Lx_02_23_for_optim = np.array(table_eRASS1_filtered_diff['L500'])#*1e42)

result_2= minimize(cost_fun_diff, x0=[1.6, 1/2.5])
B_LM_2, coeff_2 = result_2.x
print(f'---Optim result diff {result_2}')
print(f'Optimized B_LM 2 = {B_LM_2}')
print(f'Optimized Coeff 2 = {coeff_2}')

L_Maughan_2 = (uf.L_from_M_Maughan_test_B_LM(M500_eRASS1_for_optim, z_eRASS1_for_optim, B_LM=B_LM_2)*(coeff_2)) / 1e42
ratio_L_Maughan_Lx_02_23_2 = L_Maughan_2/Lx_02_23_for_optim

print(f'Mean ratio L_Maughan / Lx_02_23 = {np.mean(ratio_L_Maughan_Lx_02_23_2)}')
print(f'Standard deviation ratio L_Maughan / Lx_02_23 = {np.std(ratio_L_Maughan_Lx_02_23_2)}')

#plt.scatter(Lx_02_23_for_optim,(L_Maughan-Lx_02_23_for_optim)/(0.1*L_Maughan),alpha=0.1)
plt.figure(figsize=(8, 6))
plt.hist(ratio_L_Maughan_Lx_02_23_2 - 1, bins=50, range=[-1,1])
plt.title('Histogram of L_Maughan / Lx_02_23 pour optim 2')
plt.legend()
#plt.show()



# Troisème optimisation : cost_fun pour table low z
M500_eRASS1_for_optim = M500_eRASS1_for_optim_low_z
z_eRASS1_for_optim = z_eRASS1_for_optim_low_z
Lx_02_23_for_optim = Lx_02_23_for_optim_low_z


result_low_z = minimize(cost_fun, x0=[1.6, 1/2.5])
B_LM_optim_low_z, coeff_optim_low_z = result_low_z.x

print(f'---Optim ratio low z result {result_low_z}')
print(f'Optimized B_LM = {B_LM_optim_low_z}')
print(f'Optimized Coeff = {coeff_optim_low_z}')

L_Maughan_low_z = (uf.L_from_M_Maughan_test_B_LM(M500_eRASS1_for_optim, z_eRASS1_for_optim, B_LM=B_LM_optim_low_z)*(coeff_optim_low_z)) / 1e42
ratio_L_Maughan_Lx_02_23_low_z = L_Maughan_low_z/Lx_02_23_for_optim

print(f'Mean ratio L_Maughan / Lx_02_23 = {np.mean(ratio_L_Maughan_Lx_02_23_low_z)}')
print(f'Standard deviation ratio L_Maughan / Lx_02_23 = {np.std(ratio_L_Maughan_Lx_02_23_low_z)}')

plt.figure(figsize=(8, 6))
plt.hist(ratio_L_Maughan_Lx_02_23_low_z - 1, bins=50, range=[-1,1])
plt.title('Histogram of L_Maughan / Lx_02_23 pour optim ratio low z')
plt.legend()
#plt.show()



# Quatrième optimisation : cost_fun pour table high z
M500_eRASS1_for_optim = M500_eRASS1_for_optim_high_z
z_eRASS1_for_optim = z_eRASS1_for_optim_high_z
Lx_02_23_for_optim = Lx_02_23_for_optim_high_z

result_high_z = minimize(cost_fun, x0=[1.6, 1/2.5])
B_LM_optim_high_z, coeff_optim_high_z = result_high_z.x

print(f'---Optim ratio high z result {result_high_z}')
print(f'Optimized B_LM = {B_LM_optim_high_z}')
print(f'Optimized Coeff = {coeff_optim_high_z}')

L_Maughan_high_z = (uf.L_from_M_Maughan_test_B_LM(M500_eRASS1_for_optim, z_eRASS1_for_optim, B_LM=B_LM_optim_high_z)*(coeff_optim_high_z)) / 1e42
ratio_L_Maughan_Lx_02_23_high_z = L_Maughan_high_z/Lx_02_23_for_optim

print(f'Mean ratio L_Maughan / Lx_02_23 = {np.mean(ratio_L_Maughan_Lx_02_23_high_z)}')
print(f'Standard deviation ratio L_Maughan / Lx_02_23 = {np.std(ratio_L_Maughan_Lx_02_23_high_z)}')

plt.figure(figsize=(8, 6))
plt.hist(ratio_L_Maughan_Lx_02_23_high_z - 1, bins=50, range=[-1,1])
plt.title('Histogram of L_Maughan / Lx_02_23 pour optim ratio high z')
plt.legend()
plt.show()



'''L_Maughan_filtered = uf.L_from_M_Maughan_test_B_LM(table_eRASS1_filtered['M500']*1e13, table_eRASS1_filtered['zBest'], B_LM=B_LM_optim_2) * coeff_optim_2
ratio_L_Maughan_Lx_02_23_filtered_2 = L_Maughan_filtered/(table_eRASS1_filtered['L500'] * 1e42)

print(f'Mean ratio L_Maughan / Lx_02_23 = {np.mean(ratio_L_Maughan_Lx_02_23_filtered_2)}')
print(f'Standard deviation ratio L_Maughan / Lx_02_23 = {np.std(ratio_L_Maughan_Lx_02_23_filtered_2)}')'''


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

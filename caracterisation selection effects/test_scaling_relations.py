import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import pickle
import os
import importlib
import sys
import pandas as pd
from astropy.cosmology import Planck18


import halo_and_cluster_MF as hcmf
import useful_functions as uf

importlib.reload(hcmf) ; importlib.reload(uf)

cosmo = Planck18

# Le but de ce script est de tester les relations d'échelles fournies dans les différents papiers associés aux relevés utilisés, notamment pour vérifier si elles permettent de retrouver les mêmes valeurs

#----------------------------------------------------------------------------------------------------------------
# Test pour Gao et al. 2019 : SCUSS/SDSS/UNwise

'''

here = os.getcwd()
print('here :', here)

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

print("\nOutlier indexes for L1Mpc ratio_log10 (ratio_log10 > 1.1 or < 0.9):")
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
print(outliers_RL_null)

print('percentage of outliers L1Mpc 10% :', len(outliers_L1Mpc_10percent)/len(indexes_SCUSS)*100, '%')
print('percentage of outliers RL 10% :', len(outliers_RL_10percent)/len(indexes_SCUSS)*100, '%')
print('percentage of outliers L1Mpc extreme :', len(outliers_L1Mpc_extreme)/len(indexes_SCUSS)*100, '%')
print('percentage of outliers RL extreme :', len(outliers_RL_extreme)/len(indexes_SCUSS)*100, '%')
print('percentage of outliers L1Mpc null :', len(outliers_L1Mpc_null)/len(indexes_SCUSS)*100, '%')
print('percentage of outliers RL null :', len(outliers_RL_null)/len(indexes_SCUSS)*100, '%')

# Display the DataFrame in the terminal using pandas' to_string for better formatting
pd.set_option('display.max_rows', None)  # Show all rows (be careful with very large DataFrames)
pd.set_option('display.max_columns', None)  # Show all columns
pd.set_option('display.width', None)  # No line wrapping
print("table_SCUSS_outliers_RL_null :", table_SCUSS_outliers_RL_null.to_string(index=True))

print('')



# log10
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# First subplot: RL
axes[0].scatter(indexes_SCUSS, ratio_log10_from_RL, label='Data / retrieved', alpha=0.5)
#axes[0].scatter(indexes_SCUSS, ratio_log10_from_RL_min, label='(Data / retrieved) min', alpha=0.5)
#axes[0].scatter(indexes_SCUSS, ratio_log10_from_RL_max, label='(Data / retrieved) max', alpha=0.5)
# axes[0].axhline(y=1, color='red', label='1')
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
plt.show()

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

'''

#----------------------------------------------------------------------------------------------------------------
# Test pour eRASS1 : 
# 1er test : On utilise la relation d'échelle Lx-Mgas fournie par Bahar et al. 2022 et on compare aux données du catalogue de bulbul et al. 2024

path_eRASS1 = './DESI_eRass/eRass_Bulbul_2024_all.fit' #la première table est la table eRosita de départ, enrichie avec des informations d'association avec SDSS-Wen
table_eRASS1 = uf.fits_to_dataframe(path_eRASS1)

Mgas500_eRASS1 = np.array(table_eRASS1['Mgas500']*1e11)  ; mean_Mgas500_eRASS1 = Mgas500_eRASS1.mean()
M500_eRASS1 = np.array(table_eRASS1['M500']*1e13)  ; mean_M500_eRASS1 =  M500_eRASS1.mean()

ratio_Mgas_500_M500_eRASS1 = mean_Mgas500_eRASS1/mean_M500_eRASS1 # On va utliser ce ratio pour passer de Mgas à M500 → évidemment c'est IMPRECIS : pour la suite essayer de changer de méthode ou au moins d'estimer les erreurs
#print(ratio_Mgas_500_M500_eRASS1)



def Lx_from_Mgas500_Bahar2022(M,z):
    # z : redshift
    # Lx → erg/s
    A = 0.88 ; Aminus = 0.92-0.02 ; Aplus = 0.92+0.02
    B = 1.07 ; Bminus = 1.07-0.02 ; Bplus = 1.07+0.02
    C = 2 # Impose the gravity-only prediction for evolution
    Lx_piv = 3.20*1e43 # erg/s
    Mgas500_eRASS1_piv = 1.04*1e13 # Msun
    zpiv = 0.33 # redshift pivot value for the scaling relation

    Lx = A * (Lx_piv) * (M/Mgas500_eRASS1_piv)**(B) * (cosmo.efunc(z)/cosmo.efunc(zpiv))**(C) # cosmo.efunc(z) = E(z) = H(z)/H0

    return Lx

z_eRASS1 = np.array(table_eRASS1['zBest']) # redshift
Lx_scaled = Lx_from_Mgas500_Bahar2022(Mgas500_eRASS1, z_eRASS1) # Lx in erg/s

Lx = np.array(table_eRASS1['L500T']*1e42) # Lx in erg/s : 0.5-2 keV band
#print(Lx_scaled.mean(), Lx.mean())

# Ratio Lx_scaled / Lx
ratio_Lx_scaled_Lx = Lx_scaled/Lx


'''plt.figure(figsize=(8, 6))
plt.scatter(np.arange(len(ratio_Lx_scaled_Lx)), ratio_Lx_scaled_Lx, label='Lx_scaled / Lx', alpha=0.5)
plt.axhline(y=1, color='red', label='1')
plt.axhline(y=1.05, color='orange', label='1.05')
plt.axhline(y=0.95, color='green', label='0.95')
plt.xlabel('Index')
plt.ylabel('Lx_scaled / Lx')
plt.legend()
plt.title('Lx_scaled to Lx ratio')
plt.show()
'''
#---------------------
# 2nd test : On utilise la relation ln(Lx) - ln(m500) fournie dans chiu et al. 2022 (eq 67)

ln_Lx = np.log(Lx_scaled)

def lnLx_from_lnM500_Chiu2022(M500, z):
    A = 3.36 ; Aminus = A-0.49 ; Aplus = A + 0.53
    B = 1.44 ; Bminus = B-0.13 ; Bplus = B+0.14
    delta = -0.07 ; deltaminus = delta-0.79 ; deltaplus = delta+1.26
    gamma = -0.51 ; gammaminus = gamma-0.75 ; gammaplus = gamma+0.93
    sigma = 0.120 ; sigmaminus = sigma-0.060 ; sigmaplus = sigma+0.138
    rho_WL = 0.24 ; rho_WLminus = rho_WL-0.67 ; rho_WLplus = rho_WL+0.38
    Mpiv = 1.4 * 1e14 ; zpiv = 0.35

    result =  np.log(A) + np.log(1e43) + [ B + delta * np.log((1+z)/(1+zpiv))] * np.log(M500/Mpiv) + 2 * np.log(cosmo.efunc(z)/cosmo.efunc(zpiv)) + gamma * np.log((1+z)/(1+zpiv))

    return result


lnLxscaled = lnLx_from_lnM500_Chiu2022(M500_eRASS1, z_eRASS1)

diff_lnLxscaled_lnLx = lnLxscaled - ln_Lx

print('ratio_lnLxscaled_lnLx',np.shape(diff_lnLxscaled_lnLx))
print('lnLxscaled',np.shape(lnLxscaled))
print('ln_Lx',np.shape(ln_Lx))
print('len zeRASS1',np.shape(z_eRASS1))

print(diff_lnLxscaled_lnLx[0])

'''plt.figure(figsize=(8, 6))
plt.scatter(np.arange(len(diff_lnLxscaled_lnLx[0])), diff_lnLxscaled_lnLx[0], label='lnLxscaled - ln_Lx', alpha=0.5)
plt.axhline(y=0, color='red', label='0')
plt.axhline(y=np.log(1.05), color='orange', label='ln(1.05)')
plt.axhline(y=np.log(0.95), color='green', label='ln(1.05)')
plt.xlabel('Index')
plt.ylabel('lnLxscaled / ln_Lx')
plt.legend()
plt.title('ln_Lx_scaled to ln_Lx ratio')
plt.show()'''

#---------------------
#3 eme test : On teste la relation L-M de Maughan et al. 2018

# On calcule un première fois L_Maughan
coeff = 1/2.57 # Pour passer de bolométrique à [0.2;2.3] kev (afin de pouvoir comparer à Lx de eRASS1)
L_Maughan_ = uf.L_from_M_Maughan(M500_eRASS1,z_eRASS1) * coeff # 
Lx_02_23_ = np.array(table_eRASS1['L500']*1e42) # Lx in erg/s : 0.2-2.3keV band
ratio_L_Maughan_Lx_02_23_ = L_Maughan_/Lx_02_23_

# On filtre les outliers
table_eRASS1['ratio_L_Maughan_Lx_02_23'] = ratio_L_Maughan_Lx_02_23_
table_eRASS1_filtered = table_eRASS1[table_eRASS1['ratio_L_Maughan_Lx_02_23'] < 100]

ratio_L_Maughan_Lx_02_23 = table_eRASS1_filtered['ratio_L_Maughan_Lx_02_23']
M500_eRASS1 = np.array(table_eRASS1_filtered['M500']*1e13)  ; mean_M500_eRASS1 =  M500_eRASS1.mean()
z_eRASS1 = np.array(table_eRASS1_filtered['zBest']) # redshift

# On récupère nos valeurs 
coeff = 1/2.57
L_Maughan = uf.L_from_M_Maughan(M500_eRASS1,z_eRASS1) * coeff
Lx_02_23 = np.array(table_eRASS1_filtered['L500']*1e42) # Lx in erg/s : 0.2-2.3keV band
ratio_L_Maughan_Lx_02_23 = L_Maughan/Lx_02_23
table_eRASS1_filtered['ratio_L_Maughan_Lx_02_23'] = ratio_L_Maughan_Lx_02_23
print('mean ratio_L_Maughan_Lx_02_23 : ',np.mean(ratio_L_Maughan_Lx_02_23))
print('std ratio_L_Maughan_Lx_02_23 : ',np.std(ratio_L_Maughan_Lx_02_23))


CR_02_23 = np.array(table_eRASS1_filtered['CR500']) # CR in counts/s
CR_from_L_Maughan = uf.CR_from_LX(L_Maughan, z_eRASS1) # CR in counts/s
ratio_CR_Maughan_CR_02_23 = CR_from_L_Maughan/CR_02_23 ; ratio_CR_Maughan_CR_02_23 = ratio_CR_Maughan_CR_02_23[ratio_CR_Maughan_CR_02_23<100]
ratio_CR_outliers_0 = ratio_CR_Maughan_CR_02_23[ratio_CR_Maughan_CR_02_23 < 0.1]
print("Proportion de ratio d'éléments avec ratio CR Computed / CR Observed : ",(len(ratio_CR_outliers_0)/len(ratio_CR_Maughan_CR_02_23)))

# Plot comparison between Lx_02_23 and L_Maughan
plt.figure(figsize=(8, 6))
plt.scatter(Lx_02_23, L_Maughan, alpha=0.5, label='Data')
plt.plot([Lx_02_23.min(), Lx_02_23.max()], [Lx_02_23.min(), Lx_02_23.max()], 'r--', label='y=x')
plt.xlabel('Lx_02_23 (erg/s)')
plt.ylabel('L_Maughan (erg/s)')
plt.title('Comparison of Lx_02_23 and L_Maughan')
plt.legend()
plt.tight_layout()

'''
plt.figure(figsize=(8, 6))
plt.scatter(np.arange(len(Lx_02_23)),Lx_02_23, alpha=0.5, label='Data eRASS1')
plt.scatter(np.arange(len(Lx_02_23)),L_Maughan, alpha=0.5, label='Theory PICACS')
plt.xlabel('indexes')
plt.ylabel('L (erg/s)')
plt.title('Comparison of Lx between 0.2 and 2.3 kev and L_Maughan')
plt.legend()
plt.tight_layout()
'''
print(len(ratio_L_Maughan_Lx_02_23))

plt.figure(figsize=(8, 6))
plt.scatter(Lx_02_23,ratio_L_Maughan_Lx_02_23, alpha=0.5, label='L_Maughan/Lx_02_23')
plt.axhline(y=1, color='red')
plt.axhline(y=1.05, color='orange')
plt.axhline(y=0.95, color='green')
plt.xlabel('Lx_02_23 (erg/s)')
plt.ylabel('L (erg/s)')
plt.title('Comparison of Lx between 0.2 and 2.3 kev and L_Maughan')
plt.xscale('log')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

# Plot comparison between CR_02_23 and CR_from_L_Maughan
plt.figure(figsize=(8, 6))
plt.scatter(CR_02_23, CR_from_L_Maughan, alpha=0.5, label='Data')
plt.plot([CR_02_23.min(), CR_02_23.max()], [CR_02_23.min(), CR_02_23.max()], 'r--', label='y=x')
plt.xlabel('CR_02_23 (counts/s)')
plt.ylabel('CR_from_L_Maughan (counts/s)')
plt.title('Comparison of CR_02_23 and CR_from_L_Maughan')
plt.legend()
plt.tight_layout()

plt.figure(figsize=(12, 8))
plt.scatter(np.arange(len(ratio_CR_Maughan_CR_02_23)),ratio_CR_Maughan_CR_02_23, alpha=0.5, label='CR_Maughan/CR_02_23')
plt.axhline(y=1, color='red')
plt.axhline(y=1.05, color='orange')
plt.axhline(y=0.95, color='green')
plt.xlabel('indexes')
plt.ylabel('CR scaled / CR data (counts/s)')
plt.title('Comparison of CR between 0.2 and 2.3 kev and CR_Maughan')
plt.grid()
plt.legend()
plt.tight_layout()

plt.show()


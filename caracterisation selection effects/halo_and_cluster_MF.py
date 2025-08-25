# On définit ici les fonctions de masse des halos de matière noire et des amas de galaxie
import os
import numpy as np
import scipy
from scipy.integrate import quad
import math
import matplotlib.pyplot as plt
import importlib

from astropy.cosmology import Planck18


# Set colossus data directory before any colossus imports
#os.environ['COLOSSUS_DATA_DIR'] = r'C:\Users\ED282972\AppData\Local\colossus'''
# "COLOSSUS, a public, open-source python package for calculations related to cosmology, the large-scale structure (LSS) of matter in the universe, and the properties of dark matter halos." (Diemer 2018)
from colossus.cosmology import cosmology
from colossus.lss import mass_function

import useful_functions as uf ; importlib.reload(uf) 



#------------------------------------------------------------------------------------------------------------------
# On commence par la fonction de masse des halos de matière noire (HMF) de Murray, Robotham et Power (MRP) - Murray, Robotham, Power 2013 / Murray et al. 2021
def HMF_MRP(log10_Mmedian_Msun,log10_phimedian,alpha,beta):
    # Halo Mass Function (HMF) de Murray, Robotham et Power (MRP)

    # log10_Mmedian_Msun : log10(M_star/M_sun)
    # log10_phimedian : log10(phi_star) → Mpc^-3 dex^-1
    # alpha, beta : les 2 autres paramètres pour fitter

    def HMF_MRP_result(log10_M_Msun) :
        # log10_M_Msun : log10(M/M_sun)
        log10_M_Msun = np.array(log10_M_Msun) # Assure que c'est un tableau numpy
        M_Msun = 10**(log10_M_Msun)
        M_Msun = M_Msun #* (1/0.435) # On passe de M500 à M200 (conversion approximative voir note Reza clus_nfw.pdf) car la CMF MRP est calculée pour M200 (voir section 2.3.3 de Driver et al. 2021)
        
        phimedian = 10**(log10_phimedian)
        Mmedian_Msun = 10**(log10_Mmedian_Msun)

        result = np.log(10) * phimedian * beta * ((M_Msun * (Mmedian_Msun)**(-1))**(alpha+1)) * np.exp(-((M_Msun * (Mmedian_Msun)**(-1))**(beta)))
        return result
    
    return HMF_MRP_result


def CMF_MRP(log10_Mmedian_Msun,log10_phimedian,alpha,beta):
    # Cluster Mass Function tirée la HMF MRP (voir fonction HMF_MRP)

    def CMF_MRP_result(log10_M_Msun) :
        # log10_M_Msun : log10(M/M_sun)

        HMF = HMF_MRP(log10_Mmedian_Msun,log10_phimedian,alpha,beta)
        result = HMF(log10_M_Msun) * (1/0.85) # On considère que la masse d'un amas de galaxies provient à 85% de la matière noire (le reste étant 10% de gaz ICM et 5% de galaxies)
        
        return result
    
    return CMF_MRP_result

#---------
# Test pour ratio_M200_M500
def HMF_MRP_test(log10_Mmedian_Msun,log10_phimedian,alpha,beta):
    # Halo Mass Function (HMF) de Murray, Robotham et Power (MRP)

    # log10_Mmedian_Msun : log10(M_star/M_sun)
    # log10_phimedian : log10(phi_star) → Mpc^-3 dex^-1
    # alpha, beta : les 2 autres paramètres pour fitter

    def HMF_MRP_result_test(log10_M_Msun,ratio_M200_M500) :
        # log10_M_Msun : log10(M/M_sun)
        log10_M_Msun = np.array(log10_M_Msun) # Assure que c'est un tableau numpy
        M_Msun = 10**(log10_M_Msun)
        M_Msun = M_Msun * (ratio_M200_M500) # On passe de M500 à M200 (conversion approximative : voir note Reza clus_nfw.pdf)
        
        phimedian = 10**(log10_phimedian)
        Mmedian_Msun = 10**(log10_Mmedian_Msun)

        result = np.log(10) * phimedian * beta * ((M_Msun * (Mmedian_Msun)**(-1))**(alpha+1)) * np.exp(-((M_Msun * (Mmedian_Msun)**(-1))**(beta)))
        return result
    
    return HMF_MRP_result_test


def CMF_MRP_test(log10_Mmedian_Msun,log10_phimedian,alpha,beta):
    # Cluster Mass Function tirée la HMF MRP (voir fonction HMF_MRP)

    def CMF_MRP_result_test(log10_M_Msun,ratio_M200_M500) :
        # log10_M_Msun : log10(M/M_sun)

        HMF = HMF_MRP_test(log10_Mmedian_Msun,log10_phimedian,alpha,beta)
        result = HMF(log10_M_Msun,ratio_M200_M500) * (1/0.85) # On considère que la masse d'un amas de galaxies provient à 85% de la matière noire (le reste étant 10% de gaz ICM et 5% de galaxies)
        
        return result
    
    return CMF_MRP_result_test
#----------

def CMF_MRP_2(log10_M_Msun): 
    # L'intérêt de cette fonction est de pouvoir optimiser la valeur de fonction de masse des amas pour une masse fixée et donc en fonction des paramètres fittés
    def CMF_MRP_2_result(params): # params est un tableau de 4 paramètres : log10_Mmedian_Msun, log10_phimedian, alpha, beta
        M_Msun = 10**(log10_M_Msun)


        log10_Mmedian_Msun = params[0]; log10_phimedian = params[1] ; alpha = params[2] ; beta = params[3]
        log10_Mmedian_Msun = np.array(log10_Mmedian_Msun)  ; Mmedian_Msun = 10**(log10_Mmedian_Msun)
        log10_phimedian = np.array(log10_phimedian) ; phimedian = 10**(log10_phimedian)
        alpha = np.array(alpha)
        beta = np.array(beta)

        result = np.log(10) * phimedian * beta * ((M_Msun * (Mmedian_Msun)**(-1))**(alpha+1)) * np.exp(-((M_Msun * (Mmedian_Msun)**(-1))**(beta))) * (1/0.85) # * 1/0.85 car on considère que la DM représente 85% de la masse d'un cluster
        return result
    return CMF_MRP_2_result

def CMF_MRP_2_opposite(log10_M_Msun): 
    # Même fonction que CMF_MRP_2 mais avec résultat opposé pour maximiser en utilisant la fonction minimize de scipy.optimize
    def CMF_MRP_2_opposite_result(params): # params est un tableau de 4 paramètres : log10_Mmedian_Msun, log10_phimedian, alpha, beta
        M_Msun = 10**(log10_M_Msun)


        log10_Mmedian_Msun = params[0]; log10_phimedian = params[1] ; alpha = params[2] ; beta = params[3]
        log10_Mmedian_Msun = np.array(log10_Mmedian_Msun)  ; Mmedian_Msun = 10**(log10_Mmedian_Msun)
        log10_phimedian = np.array(log10_phimedian) ; phimedian = 10**(log10_phimedian)
        alpha = np.array(alpha)
        beta = np.array(beta)

        result = (-1) * np.log(10) * phimedian * beta * ((M_Msun * (Mmedian_Msun)**(-1))**(alpha+1)) * np.exp(-((M_Msun * (Mmedian_Msun)**(-1))**(beta))) * (1/0.85) # * 1/0.85 car on considère que la DM représente 85% de la masse d'un cluster
        return result
    return CMF_MRP_2_opposite_result



def CRF_MRP(params):  #params est un tableau de 4 paramètres : log10_Mmedian_Msun, log10_phimedian, alpha, beta
    # CRF pour "Cluster Richness Fonction" - Fonction de richesse des amas de galaxies tirée de la CMF MRP 
    # → ATTENTION : la relation M500 - Richesse est tirée de Wen et Han. 2015 (utiliseé dans Gao et al. 2020)
    def CRF_MRP_result(RL):
        log10_M_Msun = (1.08 * np.log10(np.array(RL)) - 1.37)  # relation (17) de Wen et Han. 2015 : en réalité la formule est  (1.08 +/- 0.02) * np.log10(RL) - (1.37 +/- 0.02) → essayer de prendre en compte cette dispersion pour la suite 
        log10_M_Msun += 14  # On ajoute 14 pour avoir M en unité de M_sun
        M_Msun = 10**(log10_M_Msun)

        log10_Mmedian_Msun = params[0]; log10_phimedian = params[1] ; alpha = params[2] ; beta = params[3]
        log10_Mmedian_Msun = np.array(log10_Mmedian_Msun)  ; Mmedian_Msun = 10**(log10_Mmedian_Msun)
        log10_phimedian = np.array(log10_phimedian) ; phimedian = 10**(log10_phimedian)
        alpha = np.array(alpha)
        beta = np.array(beta)

        CMF = np.log(10) * phimedian * beta * ((M_Msun * (Mmedian_Msun)**(-1))**(alpha+1)) * np.exp(-((M_Msun * (Mmedian_Msun)**(-1))**(beta))) * (1/0.85) # * 1/0.85 car on considère que la DM représente 85% de la masse d'un cluster
        coeff_multipli_list = [1.08,1.08-0.02,1.08+0.02] # plusieurs coefficients pour prendre en compte la dispersion de la relation Richesse - M500 (voir Wen et Han 2015)
        result = CMF * coeff_multipli_list[0] * (1/(np.log(10)*RL)) # calcul de dn/dRL
        result_min = CMF * coeff_multipli_list[1] * (1/(np.log(10)*RL))
        result_max = CMF * coeff_multipli_list[2] * (1/(np.log(10)*RL))

        return [result,result_min,result_max]
    return CRF_MRP_result


# !!! Pas sûr que l'impléméntation ci-dessous soit correcte (quand on trace pour SCUSS, on a pas une apparence logique et on est toujours très loin des données)
def CLF_MRP(params,z): #params est un tableau de 4 paramètres : log10_Mmedian_Msun, log10_phimedian, alpha, beta ; z est le redshift
    # CLF pour "Cluster Luminosity Fonction" - Fonction de luminosité des amas de galaxies tirée de la CMF MRP
    def CLF_MRP_result(L1Mpc):
         # Retourne un tableau de forme  len(z) * len(L1Mpc) : chaque ligne correspond à un redshift et chaque colonne à une luminosité apparente 

        L1Mpc = np.array(L1Mpc)[np.newaxis,:] # Reshape L1Mpc to a row vector
        z_ = np.array(z)[:,np.newaxis]  # Reshape z to a column vector

        RL = 0.69 * (L1Mpc**(1.32)) * (1+z_)**(2.91)  # relation (3) de Gao et al. 2019

        log10_M_Msun = 1.08 * np.log10(np.array(RL)) - 1.37  # relation (17) de Wen et Han. 2015 : en réalité la formule est  (1.08 +/- 0.02) * np.log10(RL) - (1.37 +/- 0.02) → essayer de prendre en compte cette dispersion pour la suite 
        log10_M_Msun += 14  # On ajoute 14 pour avoir M en unité de M_sun
        M_Msun = 10**(log10_M_Msun)

        log10_Mmedian_Msun = params[0]; log10_phimedian = params[1] ; alpha = params[2] ; beta = params[3]
        log10_Mmedian_Msun = np.array(log10_Mmedian_Msun)  ; Mmedian_Msun = 10**(log10_Mmedian_Msun)
        log10_phimedian = np.array(log10_phimedian) ; phimedian = 10**(log10_phimedian)
        alpha = np.array(alpha)
        beta = np.array(beta)

        CMF = np.log(10) * phimedian * beta * ((M_Msun * (Mmedian_Msun)**(-1))**(alpha+1)) * np.exp(-((M_Msun * (Mmedian_Msun)**(-1))**(beta))) * (1/0.85) # * 1/0.85 car on considère que la DM représente 85% de la masse d'un cluster
        coeff_multipli_list = [1.08,1.08-0.02,1.08+0.02] # plusieurs coefficients pour prendre en compte la dispersion de la relation Richesse - M500 (voir Wen et Han 2015)
        result = CMF * coeff_multipli_list[0] * (1/(np.log(10)*RL)) * ((0.69 * (1+z_)**2.91)*(1.32 * L1Mpc**(1.32-1))) # calcul de dn/dL # AVANR c'étatit 1.32 +1 !!! → cette correction semble avoir amélioré la comparaison avec les données théorie/données
        result_min = CMF * coeff_multipli_list[1] * (1/(np.log(10)*RL)) * ((0.69 * (1+z_)**2.91)*(1.32 * L1Mpc**(1.32-1)))
        result_max = CMF * coeff_multipli_list[2] * (1/(np.log(10)*RL)) * ((0.69 * (1+z_)**2.91)*(1.32 * L1Mpc**(1.32-1)))

        return [result,result_min,result_max]  
    return CLF_MRP_result




def CLF_log10_MRP(params,z): #params est un tableau de 4 paramètres : log10_Mmedian_Msun, log10_phimedian, alpha, beta ; z est le redshift
    # CLF pour "Cluster Luminosity Fonction" - Fonction de luminosité des amas de galaxies tirée de la CMF MRP
    # On mène les calculs à z fixés !
    def CLF_log10_MRP_result(log10_L1Mpc):
         # Retourne un tableau de forme  len(z) * len(log10_L1Mpc) : chaque ligne correspond à un redshift et chaque colonne à une luminosité  

        log10_L1Mpc = np.array(log10_L1Mpc)[np.newaxis,:] # Reshape L1Mpc to a row vector
        z_ = np.array(z)[:,np.newaxis]  # Reshape z to a column vector

        L1Mpc = 10**(log10_L1Mpc) # We will need that quantity for the computation of the CLF using the CMF
        RL = 0.69 * (L1Mpc**(1.32)) * (1+z_)**(2.91)  # Idem (from relation (3) de Gao et al. 2019)

        log10_M_Msun = 1.08 * (np.log10(0.69) + 1.32 * log10_L1Mpc + 2.91 * np.log10(1+z_)) - 1.37 
        # On injecte la relation (3) de Gao et al. 2019 dans la relation (17) de Wen et Han. 2015 : en réalité la formule est  (1.08 +/- 0.02) * np.log10(RL) - (1.37 +/- 0.02) → essayer de prendre en compte cette dispersion pour la suite 
        log10_M_Msun += 14  # On ajoute 14 pour avoir M en unité de M_sun
        M_Msun = 10**(log10_M_Msun)

        log10_Mmedian_Msun = params[0]; log10_phimedian = params[1] ; alpha = params[2] ; beta = params[3]
        log10_Mmedian_Msun = np.array(log10_Mmedian_Msun)  ; Mmedian_Msun = 10**(log10_Mmedian_Msun)
        log10_phimedian = np.array(log10_phimedian) ; phimedian = 10**(log10_phimedian)
        alpha = np.array(alpha)
        beta = np.array(beta)

        CMF = np.log(10) * phimedian * beta * ((M_Msun * (Mmedian_Msun)**(-1))**(alpha+1)) * np.exp(-((M_Msun * (Mmedian_Msun)**(-1))**(beta))) * (1/0.85) # * 1/0.85 car on considère que la DM représente 85% de la masse d'un cluster
        coeff_multipli_list = [1.08,1.08-0.02,1.08+0.02] # plusieurs coefficients pour prendre en compte la dispersion de la relation Richesse - M500 (voir Wen et Han 2015)
        result = CMF * coeff_multipli_list[0] * (1/(np.log(10)*RL)) * ((0.69 * (1+z_)**2.91)*(1.32 * L1Mpc**(1.32-1))) * (np.log(10)*L1Mpc) # calcul de dn/dL
        result_min = CMF * coeff_multipli_list[1] * (1/(np.log(10)*RL)) * ((0.69 * (1+z_)**2.91)*(1.32 * L1Mpc**(1.32-1))) * (np.log(10)*L1Mpc)  # calcul de dn/dL
        result_max = CMF * coeff_multipli_list[2] * (1/(np.log(10)*RL)) * ((0.69 * (1+z_)**2.91)*(1.32 * L1Mpc**(1.32-1))) * (np.log(10)*L1Mpc)  # calcul de dn/dL
        return [result,result_min,result_max]  
    return CLF_log10_MRP_result


def CLsunF_log10_MRP(params,z): 
    # Il s'agit de la 'Cluster Luminosity Function' pour une unité de luminosité en Lsun
    # params est un tableau de 4 paramètres : log10_Mmedian_Msun, log10_phimedian, alpha, beta ; z est le redshift
    # On mène les calculs à z fixés !
    def CLsunF_log10_MRP_result(log10_Lsun_L1Mpc):
        #Retourne un tableau de forme  len(z) * len(log10_Lsun_L1Mpc) : chaque ligne correspond à un redshift et chaque colonne à une luminosité

        log10_Lsun_L1Mpc = np.array(log10_Lsun_L1Mpc)[np.newaxis,:] # Reshape L1Mpc to a row vector : log10(L1Mpc/Lsun)
        z_ = np.array(z)[:,np.newaxis]  # Reshape z to a column vector

        log10_L1Mpc = log10_Lsun_L1Mpc - np.log10((1.08 * 10**10)*10**(0.4*1.16*z_)) # ← log10_L1Mpc = log10(L/L*) ; avec L*(z=0) = 1.08 * 10**10 * Lsun
        # ATTENION ! La correction écrase tout pour des z trop grands

        log10_M_Msun = 1.08 * (np.log10(0.69) + 1.32 * log10_L1Mpc + 2.91 * np.log10(1+z_)) - 1.37 
        # On injecte la relation (3) de Gao et al. 2019 dans la relation (17) de Wen et Han. 2015 : en réalité la formule est  (1.08 +/- 0.02) * np.log10(RL) - (1.37 +/- 0.02) → essayer de prendre en compte cette dispersion pour la suite 
        log10_M_Msun += 14  # On ajoute 14 pour avoir M en unité de M_sun
        M_Msun = 10**(log10_M_Msun)

        log10_Mmedian_Msun = params[0]; log10_phimedian = params[1] ; alpha = params[2] ; beta = params[3]
        log10_Mmedian_Msun = np.array(log10_Mmedian_Msun)  ; Mmedian_Msun = 10**(log10_Mmedian_Msun)
        log10_phimedian = np.array(log10_phimedian) ; phimedian = 10**(log10_phimedian)
        alpha = np.array(alpha)
        beta = np.array(beta)

        CMF = np.log(10) * phimedian * beta * ((M_Msun * (Mmedian_Msun)**(-1))**(alpha+1)) * np.exp(-((M_Msun * (Mmedian_Msun)**(-1))**(beta))) * (1/0.85) # * 1/0.85 car on considère que la DM représente 85% de la masse d'un cluster
        coeff_multipli_list = [1.08,1.08-0.02,1.08+0.02] # plusieurs coefficients pour prendre en compte la dispersion de la relation Richesse - M500 (voir Wen et Han 2015)
        result = CMF * coeff_multipli_list[0] * (1.32*log10_L1Mpc) # calcul de dn/dL
        result_min = CMF * coeff_multipli_list[1] * (1.32*log10_L1Mpc) # calcul de dn/dL
        result_max = CMF * coeff_multipli_list[2] * (1.32*log10_L1Mpc) # calcul de dn/dL

        return [result,result_min,result_max]  
    return CLsunF_log10_MRP_result
#----------------------------------------------------------------------------------------------------------------
# On définit ici la fonction de masse de halo foournie par Bocquet et al. 2016
# A partir de maintenant on va saider du package COLOSSUS


def HMF_Bocquet_2016(M,z):
    # Halo mass function
    # M : range of masses in Msun
    # z : redshift
    cosmo_ = cosmology.setCosmology('planck18', persistence='') # persistence='' ralentit le processus mais empêche des problèmes d'accès à la mémoire cache du package colossus
    result =  mass_function.massFunction(M, z, mdef = '500c', model = 'bocquet16', q_out = 'dndlnM') # dndlnM → return dn/dln(M)

    return result

def CMF_Bocquet_2016(M,z):
    # Cluster mass function 
    # M : range of masses in Msun
    # z : redshift
    cosmo_ = cosmology.setCosmology('planck18', persistence='')
    result =  mass_function.massFunction(M, z, mdef = '500c', model = 'bocquet16', q_out = 'dndlnM')* (1/0.85)  # On considère que la masse d'un amas de galaxies provient à 85% de la matière noire (le reste étant 10% de gaz ICM et 5% de galaxies)

    return result # (dn/dln(M))(M)

def CMF_Bocquet_2016_log10(M,z):
    # Cluster mass function 
    # M : range of masses in Msun
    # z : redshift
    cosmo_ = cosmology.setCosmology('planck18', persistence='')
    result =  mass_function.massFunction(M, z, mdef = '500c', model = 'bocquet16', q_out = 'dndlnM')* (1/0.85)  # On considère que la masse d'un amas de galaxies provient à 85% de la matière noire (le reste étant 10% de gaz ICM et 5% de galaxies)

    return result * np.log(10) # (dn/dlog10(M))(M)
 
def CLxF_eRASS_Bocquet_2016_DEPRECATED(Lx, z):
    # DEPRECATED
    # Cluster X-ray luminosity function adapted for eRASS1 (Bulbul et al. 2024, Ghirardini et al. 2024) and using Bocquet et al. 2016 CMF
    # Lx : range of luminosities in erg/s
    # z : redshift
    # Attention on considère un z fixe ici ! (Faire une surface en 3D pour la suite ?)

    # On commence par calculer une APPROXIMATION de M500(CR(Lx)) avec CR le Count Rate de photons entre 0.2 et 2.3 kev (Lx correspond à la même bande)
    # Lx → CR
    cosmo = Planck18
    alpha = 0.42 # alpha =/= 1 est légitime ? discussion Réza !
    DLz = (cosmo.luminosity_distance(z) * (3.086 * 1e22)) # Distance lumineuse en m (attention c'est la même distance lumineuse pour chaque calcul dû au z fixe : grosse approximation !)
    SeROSITA = 2451 * 1e-4 # Surface effective totale de collection des photons avec eROSITA en m^2 (calculé dans le notebook)
    E_peak_kev = np.round(2.821 * 0.5740, 3) #localisation du pic dans le spectre du corps noir pour la température moyenne des amas dans eRASS1 (temp moyenne = 0.5470 kev + loi de déplacement de Wien)
    E_peak = E_peak_kev * (1.6022 * 1e-9) # unit : erg

    CR_from_LX = alpha * (Lx/E_peak) * (SeROSITA/(4*np.pi*(DLz)**2)) 

    # CR → M500 (Ghirardini et al. 2024)
        # Pivot values
    zp = 0.35  # Redshift pivot
    M500p = 2*(1e14) #Masse pivot - M_sun
    CRp = 0.1 # Count rate pivot - cts/s

    Ax = 0.64 ; Ax_min = Ax - 0.06 ; Ax_max = Ax + 0.04

    Dx = -2 ; Ex = 2
    Gx = 0.29 ; Gx_min = Gx - 0.13 ; Gx_max = Gx + 0.12
    ex_z = Dx * np.log(cosmo.luminosity_distance(z)/cosmo.luminosity_distance(zp)) + Ex * np.log(cosmo.efunc(z)/cosmo.efunc(zp)) + Gx * np.log((1+z)/(1+zp))

    Bx = 1.38 ; Bx_min = Bx - 0.04 ; Bx_max = Bx + 0.03
    Fx = -0.33 ; Fx_min = Fx - 0.12 ; Fx_max = Fx + 0.12
    bx_z = Bx + Fx * np.log((1+z)/(1+zp))

    ln_M500 = np.log(M500p) + (1/bx_z) * (np.log(np.asarray(CR_from_LX/CRp)) - np.log(Ax) - ex_z) # np.asarray : transform quantitities inside into floats → avoid error with argument of log with a dimension
    M500 = np.exp(ln_M500)  # M500 en M_sun
    print('M500 = ', M500)

    #Maintenant qu'on dispose d'une approximation de M500, on peut calculer notre ClxF en partant de la CMF de Bocquet 2016
    dn_dlnM = CMF_Bocquet_2016(M500,z)

    # On calcule les éléments différentiels dlnM_dlnCR et dlnCR_dlnLx
    dlnM_dlnCR = bx_z

    dlnCR_dlnLx = 1 # voir mail reza 09/08/2025

    #On a notre résultal final !
    dn_dlnLx = dn_dlnM * dlnM_dlnCR * dlnCR_dlnLx

    return dn_dlnLx


def CLxF_eRASS_Bocquet_2016(Lx , z):
    # Cluster X-ray luminosity function adapted for eRASS1 (Maughan et al. 2013) and using Bocquet et al. 2016 CMF
    # Lx : Cluster Luminosity in erg/s, between 0.2 and 2.3 keV
    # z : redshift

    # Renvoie (dn/dlog10(Lx))(Lx(M500))

    # On commence par calculer (dn/dlog10(M500))(M500) isssu de Bocquet et al. 2016
    M500 = uf.M_from_Lx(Lx,z) #→ utilise Maughan et al. 2013
    dn_dlog10M500 = CMF_Bocquet_2016_log10(M500, z)

    # On calcule l'élément différentiel
    a = -0.495 ; b = 1.539
    B_LM = a * z + b

    dlog10M500_dlog10Lx = 1/B_LM

    # On a notre résultat final :
    dn_dlog10Lx = dn_dlog10M500 * dlog10M500_dlog10Lx

    return dn_dlog10Lx




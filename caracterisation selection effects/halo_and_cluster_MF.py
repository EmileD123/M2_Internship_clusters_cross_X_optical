# On définit ici les fonctions de masse des halos de matière noire et des amas de galaxie

import numpy as np
import scipy
from scipy.integrate import quad
import math
import matplotlib.pyplot as plt

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

def hello_world():
    print("Hello, world! version 5")







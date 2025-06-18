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



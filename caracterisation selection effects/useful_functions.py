import numpy as np
import scipy
from scipy.integrate import quad
import math
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18
from astropy.io import fits
import pandas as pd


cosmo = Planck18

# Calcule la densité volumique moyenne des amas de galaxies en intégrant la fonction de masse des amas (CMF) sur une plage de masses donnée
def CMF_to_mean_density(CMF,bornes_log10_M_Msun):
    # CMF : fonction de masse des amas de galaxies considérée
    # bornes_log10_M_Msun : tuple contenant les bornes inférieure et supérieure de l'intervalle de masse en log10(M/M_sun)

    result = quad(CMF, bornes_log10_M_Msun[0], bornes_log10_M_Msun[1])
    return result[0]  # Le résultat de quad est un tuple, le premier élément est la valeur de l'intégrale



# Calcule le volume (en MPC^3) couvert par un relevé à l'aide de sa surface angulaire et de l'intervalle de redshift
def Volume_survey(Angular_surface_deg2,bornes_z):
    # Angular_surface_deg2 : surface angulaire du relevé en degrés carrés
    # bornes_z : tuple contenant les bornes inférieure et supérieure de l'intervalle de redshift

    Angular_surface_sr = Angular_surface_deg2 * (np.pi / 180)**2  # Conversion de degrés carrés en stéradians
    D1 = cosmo.comoving_distance(bornes_z[0]).value ; D2 = cosmo.comoving_distance(bornes_z[1]).value  # On prend la distance comobile pout le calcul du volume, ce qui permet d'ignorer l'expansion de l'univers,
    #c'est ce qui est utlisée pour le calcul des HMF (voir par exemple Kosiba et al. 2024, section 2.1)
    volume = (D2**3 - D1**3) * Angular_surface_sr / 3  # Volume d'un cône tronqué en MPC^3
    return volume

# Transforme certaines extensions d'un fichier FITS en un fichier FITS ASCII
def fits_binary_to_fits_ascii_list(fits_file, liste_indexes, name_list):  # liste_indexes indique les extensions du fits à lire ; name_list correspond aux noms des fichiers de sortie
    for x, x_name in zip(liste_indexes, name_list):  # Use zip to iterate over both lists simultaneously
        with fits.open(fits_file) as hdul:
            data = hdul[x].data
            data_table = Table(data)  # Convertit le tableau binaire en un tableau astropy
            data_table.write(x_name, overwrite=True, format='ascii')  # Écrit le tableau dans un fichier FITS ASCII

#Transforme un FITS ascii en un DataFrame pandas
def ascii_to_dataframe(ascii_file, delimiter=' '):
    data_df = pd.read_csv(ascii_file, delimiter=delimiter)
    return data_df

# Transfome la première extension d'un fichier FITS en un DataFrame pandas
def fits_to_dataframe(fits_file):
    with fits.open(fits_file) as hdul:
        data = hdul[1].data 
        data_df = pd.DataFrame(data)
    return data_df




# Transforme certaines extensions d'un fichier FITS en un DataFrame pandas(→ pour que ça fonctionne, il faut un fits ascii et pas un fits binaire)
def fits_to_dataframe_list(fits_file, liste_indexes):  # liste_indexes indique les extensions du fits à lire
    df_list = []  # Initialize the list outside the loop
    with fits.open(fits_file) as hdul:
        for x in liste_indexes:
            data = hdul[x].data  # Select the data in the xth extension
            df_list.append(pd.DataFrame(data))
            print(df_list)
    return df_list  # Return the entire list of DataFrames





# Relations importantes dans le papier Gao et al. 2019 (SCUSS/SDSS/unWISE)

def L1Mpc_to_log10_M_Msun(L1Mpc, z):
    # Relation entre la luminosité L1Mpc et la M500 des amas de galaxies, en passant par le calcul de la richesse
    # L1Mpc : "is the total r-band luminosity of cluster members within a distance of 1 Mpc from the BCG" → luminosité absolue et non apparente
    # z : redshift
    # Retourne un tableau de forme  len(z) * len(L1Mpc) : chaque ligne correspond à un redshift et chaque colonne à une luminosité apparente 

    L1Mpc = np.array(L1Mpc)[np.newaxis,:] # Reshape L1Mpc to a row vector
    z = np.array(z)[:,np.newaxis]  # Reshape z to a column vector

    richness = 0.69 * (L1Mpc**1.32) * (1 + z)**2.91  # relation (3) de Gao et al. 2019

    log10_M_Msun = 1.08 * np.log10(richness) - 1.37  # relation (17) de Wen et al. 2015
    log10_M_Msun += 14  # On ajoute 14 pour avoir M en unité de M_sun

    return log10_M_Msun   # Returns an array of size len(z)*len(L1Mpc)


def RL_to_log10_M_Msun(richness):
    # !!! ATTENTION : La fonction renvoie 3 tableaux de masses (la moyenne  et les bornes inférieure et supérieure) pour prendre en compte la dispersion)
    # Relation entre la richesse et la M500 des amas de galaxies
    # Retourne la masse en log10(M/M_sun)

    log10_M_Msun = 1.08 * np.log10(richness) - 1.37  # relation (17) de Wen et al. 2015 : en réalité la formule est  (1.08 +/- 0.02) * np.log10(richness) - (1.37 +/- 0.02) → essayer de prendre en compte cette dispersion pour la suite
    # M_Msun est pour l'instant en unité de 10^14 M_sun
    log10_M_Msun += 14  # On ajoute 14 pour avoir M en unité de M_sun
    
    log10_M_Msun_min =  (1.08 - 0.02) * np.log10(richness) - (1.37+0.02) + 14
    log10_M_Msun_max =  (1.08 + 0.02) * np.log10(richness) - (1.37-0.02) + 14

    return [log10_M_Msun,log10_M_Msun_min,log10_M_Msun_max]


import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import fsolve
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
    # bornes_z : [a,b] → tuple contenant les bornes inférieure et supérieure de l'intervalle de redshift

    Angular_surface_sr = Angular_surface_deg2 * (np.pi / 180)**2  # Conversion de degrés carrés en stéradians
    D1 = cosmo.comoving_distance(bornes_z[0]).value ; D2 = cosmo.comoving_distance(bornes_z[1]).value  # On prend la distance comobile pout le calcul du volume, ce qui permet d'ignorer l'expansion de l'univers,
    #c'est ce qui est utlisée pour le calcul des HMF (voir par exemple Kosiba et al. 2024, section 2.1)
    volume = (D2**3 - D1**3) * Angular_surface_sr / 3  # Volume d'un cône tronqué en MPC^3
    return volume

# Idem mais fonctionne pour de multiples bornes de z (typiquement dans le cas où on à un zmax pour chaque bin de luminosité)
def Volume_survey_array(Angular_surface_deg2,bornes_z):
    # Angular_surface_deg2 : surface angulaire du relevé en degrés carrés
    # bornes_z : [a,[b]] tuple contenant les bornes inférieure et supérieure de l'intervalle de redshift ; [b] de taille > 1
    volume = []
    Angular_surface_sr = Angular_surface_deg2 * (np.pi / 180)**2  # Conversion de degrés carrés en stéradians
    D1 = cosmo.comoving_distance(bornes_z[0]).value 
    for zmax in bornes_z[1]:
        D2 = cosmo.comoving_distance(zmax).value 
        current_volume = (D2**3 - D1**3) * Angular_surface_sr / 3
        volume.append(current_volume)
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




#----------------------------------------------------------------------------------------------------------------
# IMPORTANT : Relations VALABLES dans le papier Gao et al. 2019 (SCUSS/SDSS/unWISE)
# (générlisables aux autres relevés ?)

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

def L1Mpc_to_log10_M_Msun_test(L1Mpc, z): # Dédié uniquement au test dans test_scaling_relations.py
    # Relation entre la luminosité L1Mpc et la M500 des amas de galaxies, en passant par le calcul de la richesse
    # L1Mpc : "is the total r-band luminosity of cluster members within a distance of 1 Mpc from the BCG" → luminosité absolue et non apparente
    # z : redshift
    # Retourne un tableau de forme  len(z) * len(L1Mpc) : chaque ligne correspond à un redshift et chaque colonne à une luminosité apparente 

    L1Mpc = np.array(L1Mpc)
    z = np.array(z)

    # If L1Mpc and z are not the same length, raise an error
    if L1Mpc.shape[0] != z.shape[0]:
        raise ValueError("L1Mpc and z must have the same length to extract the diagonal.")

    richness = 0.69 * (L1Mpc**1.32) * (1 + z)**2.91  # relation (3) de Gao et al. 2019

    log10_M_Msun = 1.08 * np.log10(richness) - 1.37  # relation (17) de Wen et al. 2015
    #log10_M_Msun += 14  # On ajoute 14 pour avoir M en unité de M_sun

    return log10_M_Msun   # Returns only the diagonal (element-wise) result

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

def RL_to_log10_M_Msun_test(richness): # Dédié uniquement au test dans test_scaling_relations.py (on ne convertit pas en M_sun en ajoutant 14)
    # !!! ATTENTION : La fonction renvoie 3 tableaux de masses (la moyenne  et les bornes inférieure et supérieure) pour prendre en compte la dispersion)
    # Relation entre la richesse et la M500 des amas de galaxies
    # Retourne la masse en log10(M/M_sun)

    log10_M_Msun = 1.08 * np.log10(richness) - 1.37  # relation (17) de Wen et al. 2015 : en réalité la formule est  (1.08 +/- 0.02) * np.log10(richness) - (1.37 +/- 0.02) → essayer de prendre en compte cette dispersion pour la suite
    # M_Msun est pour l'instant en unité de 10^14 M_sun
    
    #log10_M_Msun += 14  # On ajoute 14 pour avoir M en unité de M_sun
    
    log10_M_Msun_min =  (1.08 - 0.02) * np.log10(richness) - (1.37+0.02) #+ 14
    log10_M_Msun_max =  (1.08 + 0.02) * np.log10(richness) - (1.37-0.02) #+ 14

    return [log10_M_Msun,log10_M_Msun_min,log10_M_Msun_max]



def inverse_luminosity_distance(z_array):
    # z_array : [start,stop,num]
    # Return a function that approximates the inverse of DL(z) on the range defined by z_array
    start = z_array[0] ; stop = z_array[1] ; num = z_array[2] 
    range_z = np.linspace(start, stop, num)
    DL = cosmo.luminosity_distance(range_z).value * 1e6
    # Now, return a function that approximates the inverse of DL(z)
    # We'll use interpolation to create an approximate inverse function
    from scipy.interpolate import interp1d
    # DL is strictly increasing with z, so we can invert it
    DL_sorted = np.array(DL)
    z_sorted = np.array(range_z)
    # Create the inverse function: given DL, return z
    inverse_func = interp1d(DL_sorted, z_sorted, kind = 'linear', bounds_error=False)
    return inverse_func


def inverse_luminosity_distance_old(Dl):
    # Deprecated
    # Dl est en pc
    # Compute the inverse of the luminosity distance function to find z for a given Dl.
    Dl = np.array(Dl)
    size_array = np.size(Dl)
    def func(z):
        return (cosmo.luminosity_distance(z).value * 1e6) - Dl # en pc
    
    # Use fsolve to find the root of the function
    z_initial_guess = np.ones(size_array) * 0.5  # Initial guess for the redshift
    z_solution = fsolve(func, z_initial_guess)
    return z_solution  # Return the first element of the solution array, which is the redshift z we search for


def L_star_L_sun_for_SCUSS(table):
    # Ajoute une colonne log10(L/L_sun) à la table issu de Gao et al. 2020 → uniquement celle là car ils ont un calcul de L* spécifique venant de Blanton et al. 2003
    # 2 limitations importantes pour cette relation : 1) La formule est ajustée empiriquement dans la plage de redshift : 0.02<z<0.22
    #                                                 2) L'ensemble des galaxies est supposée évolué identiquement

    Lstar_z_0 = 1.08 * 1e10 # unité : L_sun
    Q = 1.16

    log10_L1mpc = np.log10(table['L1Mpc']) # log10(L/L*)
    z = np.array(table['z'])

    L_star = Lstar_z_0 * 10**(0.4 * Q * z) # unité : L_sun
    result = log10_L1mpc - np.log10(1/L_star)

    return result

#----------------------------------------------------------------------------------------------------------------    
# Fonctions anciennement dans le notebook (faites par CURSOR)

# Fonctions pour raccourcir le code
def compute_histogram_density(data, bins=100, range=None, volume=1.0): #Attention : on ne ressort ici 
    """
    Compute histogram, bin centers, bin width, and density for a given data array.

    Args:
        data (array-like): Input data (e.g., log-mass, richness, luminosity).
        bins (int): Number of bins.
        range (tuple): Range for the histogram.
        volume (float or array-like): Survey volume for normalization.

    Returns:
        bin_centers (np.ndarray): Centers of the bins.
        counts (np.ndarray): Counts per bin.
        bin_width (float): Width of each bin.
        density (np.ndarray): Density per bin (counts/bin_width/volume).
    """
    counts, bin_edges = np.histogram(data, bins=bins, range=range)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width = np.diff(bin_edges)[0]  # Assuming uniform bin width
    density = (counts/bin_width)/volume
    return bin_centers, counts, bin_width, density

def compute_poisson_errors(counts, bin_width, normalize=False):
    """
    Compute Poisson errors for histogram bins, optionally normalized.

    Args:
        counts (np.ndarray): Counts per bin.
        bin_width (float): Bin width.
        normalize (bool): If true, normalize errors by (total * bin_width)

    Returns:
        errors (np.ndarray): Poisson errors (normalized if total is given).
    """
    errors = np.sqrt(counts)
    if normalize:
        total = np.sum(counts)
        errors = errors / (total * bin_width)
    return errors

def process_subtables_histogram(subtables, column, bins, range, volume, label_prefix, color_map):
    """
    Compute and plot histograms/densities for a list of subtables (e.g., redshift bins).

    Args:
        subtables (list of DataFrame): List of subtables.
        column (str): Column name to process.
        bins (int): Number of bins.
        range (tuple): Histogram range.
        volume (array-like): Survey volume.
        label_prefix (str): Prefix for legend labels.
        color_map (callable): Function mapping index to color.
    """
    results = []
    for i, subtable in enumerate(subtables):
        data = subtable[column]
        bin_centers, counts, bin_width, density = compute_histogram_density(
            data, bins=bins, range=range, volume=volume[i]
        )
        norm_counts = counts / (np.sum(counts) * bin_width)
        errors = compute_poisson_errors(counts, bin_width)
        norm_errors = compute_poisson_errors(counts, bin_width, normalize=True)
        
        results.append({
            'bin_centers': bin_centers,
            'counts': counts,
            'norm_counts': norm_counts,
            'density': density,
            'errors': errors,
            'norm_errors': norm_errors,
            'color': color_map(i / (len(subtables) - 1))
        })
    return results

def add_apparent_luminosity(df, Labs_col, z_col, cosmo, new_col='Lapp'):
    """
    Add a column for apparent luminosity to a DataFrame.

    Args:
        df (pd.DataFrame): Input DataFrame.
        Labs_col (str): Name of the absolute luminosity column.
        z_col (str): Name of the redshift column.
        cosmo (astropy.cosmology): Cosmology object for distance calculation.
        new_col (str): Name for the new column.

    Returns:
        pd.DataFrame: DataFrame with new column added.
    """
    df = df.copy()
    Lapp = []
    for _, row in df.iterrows():
        Labs = row[Labs_col]
        z = row[z_col]
        DL = (cosmo.luminosity_distance(z).value)*1e6  # Distance lumineuse en pc
        Lapp.append(((10/DL)**2) * Labs)  # Conversion de L1Mpc à Lapp
    df[new_col] = Lapp
    return df

    #---------------------------------------------------------------------------------------------------

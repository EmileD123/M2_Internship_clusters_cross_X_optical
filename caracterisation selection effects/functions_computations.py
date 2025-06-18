import numpy as np
import scipy
from scipy.integrate import quad
import math
import matplotlib.pyplot as plt
from astropy.cosmology import Planck13

cosmo = Planck13

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
    D1 = cosmo.angular_diameter_distance(bornes_z[0]).value ; D2 = cosmo.angular_diameter_distance(bornes_z[1]).value  # On prend la distance angulaire correpondant aux bornes de redshift (pourquoi pas un autre type de deistance, comme la comobile ? ← comme la comobile)
    volume = (D2**3 - D1**3) * Angular_surface_sr / 3  # Volume d'un cône tronqué en MPC^3
    return volume
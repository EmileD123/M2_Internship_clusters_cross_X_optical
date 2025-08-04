import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import math
import pickle
import os
from scipy.optimize import minimize, Bounds
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from astropy.cosmology import Planck18
import importlib
from scipy.interpolate import UnivariateSpline
from scipy import integrate

import halo_and_cluster_MF as hcmf
import useful_functions as uf

# Permet de recharger le module hcmf et uf
importlib.reload(hcmf) ; importlib.reload(uf)

#import data------------------------------------------------
path_SCUSS = '../SCUSS_DESI/SCUSS_Gao_2020.fit' #la première table est la table eRosita de départ, enrichie avec des informations d'association avec SDSS-Wen
table_SCUSS = uf.fits_to_dataframe(path_SCUSS)
M500_SCUSS = np.log10(np.array(table_SCUSS['M500'])*1e13)

path_eRASS1 = '../DESI_eRass/eRass_Bulbul_2024.fit' #la première table est la table eRosita de départ, enrichie avec des informations d'association avec SDSS-Wen
table_eRASS1 = uf.fits_to_dataframe(path_eRASS1)
M500_eRASS1 = np.log10(np.array(table_eRASS1['M500'])*1e13)
#------------------------------------------------------

#DATA
## SCUSS/SDSS/unWISE
bornes_z_SCUSS = [0.05,0.65]
angular_area_SCUSS_deg2 = 3700 # Surface couverte par SCUSS/SDSS/unWISE en degrés carrés
volume_SCUSS = uf.Volume_survey(angular_area_SCUSS_deg2,bornes_z_SCUSS)

bin_centers, counts, bin_width_SCUSS, n_SCUSS = uf.compute_histogram_density(
    data=M500_SCUSS,
    bins=100,
    volume=volume_SCUSS
)
print("Bin width:", bin_width_SCUSS)


params_GSR = [14.13,-3.96,-1.68,0.63] # valeurs moyennes 
bounds_params_GSR = [(14.13-0.40,14.13+0.43),(-3.96-0.82,-3.96+0.55),(-1.68-0.24,-1.68+0.21),(0.63-0.11,0.63+0.25)]

CMF_GSR_func = hcmf.CMF_MRP_test(params_GSR[0],params_GSR[1],params_GSR[2],params_GSR[3])

log10_M_Msun = np.linspace(12.7,15.5,1000)
M_Msun = 10**log10_M_Msun
log10_M_Msun = log10_M_Msun.tolist()






#ratio_M200_M500 = 0.1

def test_ratio_M200_M500(ratio_M200_M500):
    estimation_data = np.array(n_SCUSS / CMF_GSR_func(bin_centers,ratio_M200_M500))
    proxy_area = np.sum(estimation_data)
    return proxy_area

def constraint_max_estimation_data(ratio):
    # Compute estimation_data for this ratio
    estimation_data = np.array(n_SCUSS / CMF_GSR_func(bin_centers, ratio))
    # The constraint is that all elements must be ≤ 1
    # So, for minimize, we want 1 - estimation_data ≥ 0 for all elements
    return 1 - estimation_data

def objective(ratio):
    # Negative because we want to maximize
    return -test_ratio_M200_M500(ratio[0])

# Set reasonable bounds for ratio_M200_M500, e.g., between 0.5 and 2
bounds = Bounds([1/0.9], [1/0.1])

# Initial guess
x0 = [0.15]

# Define the constraint in the form required by scipy.optimize
cons = {'type': 'ineq', 'fun': lambda x: constraint_max_estimation_data(x[0])}

result = minimize(objective, x0, bounds=bounds, constraints=cons)

print("Optimal ratio_M200_M500:", result.x[0])
print("Maximum proxy_area:", -result.fun)


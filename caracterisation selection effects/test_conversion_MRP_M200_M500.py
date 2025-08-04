import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import math
import pickle
import os
from scipy.optimize import minimize
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from astropy.cosmology import Planck18
import importlib
from scipy.interpolate import UnivariateSpline
from scipy import integrate

import halo_and_cluster_MF as hcmf
import useful_functions as uf

# Calcul classique de la fonction de masse (dn/dlog10(M200))[M200]

params_GSR = [14.13,-3.96,-1.68,0.63] # valeurs moyennes : pour z = 0.1 !
CMF_GSR_func_M200 = hcmf.CMF_MRP(params_GSR[0],params_GSR[1],params_GSR[2],params_GSR[3])
log10_M_Msun_M200 = np.linspace(12.7,15.5,1000) # Les fits desquels sont issus les paramètres sont valables pour des masses comprises entre 10^12.7 et 10^15.5 Msun
M_Msun_200 = 10**log10_M_Msun_M200
log10_M_Msun_M200 = log10_M_Msun_M200.tolist()
CMF_GSR_M200 = CMF_GSR_func_M200(log10_M_Msun_M200)

'''plt.figure(figsize=(10, 6))
plt.plot(log10_M_Msun, CMF_GSR_M200, label='CMF GSR M200', color='blue')
plt.xlabel('log10(M/Msun)')
plt.ylabel('dn/dlog10(M200) [Mpc^-3]')
plt.title('Fonction de Masse (CMF) pour M200')
plt.yscale('log')
plt.grid()
plt.legend()
plt.show()'''
 
# Calcul d'un array de log10(M500/Msun) à partir d'un array d'un array de log10(M200/Msun) thanks to Duffy et al. 2008



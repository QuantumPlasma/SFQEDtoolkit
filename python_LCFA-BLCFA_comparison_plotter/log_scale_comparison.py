import numpy as np
import matplotlib.pyplot as plt

from scipy.special import kv
from scipy.integrate import quad
import scipy.integrate as integrate

# --------------------------------------------------------- #
# ANALYTICAL 
# --------------------------------------------------------- #
# Constants (SI units)
# --------------------------------------------------------- #

colors = ['#0036A6', '#FF6319', '#00933C', '#FC0009', '#60269E', '#996633',
          '#B933AD', '#808183', '#FFBE00', '#6CBE45', '#A7A9AC', '#000000']
cnames = ['blue', 'orange', 'green', 'red', 'violet', 'brown', 'magenta',
          'grey', 'yellow', 'lightgreen', ' lightgrey']

alpha = 1./137. 
HBar = 1.0545718e-34 # (*m^2kg/s*)
m = 9.10938356e-31 # (*kg*)
c = 299792458 # (* m/s *)
e = 1.60217662e-19 # (*coulombs*)
ptive_inf = float('inf')

wr_dGamma_coeff = (alpha * m * c * c) /( np.sqrt(3) * np.pi * HBar)

def K1By3(y):
  return kv(1./3., y)

def Wr_dGamma(gamma_en, chi_e, ele_en):
  v = (2.*gamma_en)/(3.*chi_e*(ele_en-gamma_en))
  term1 = (ele_en - gamma_en)/(ele_en*ele_en*ele_en)
  term2 = ( 1. +  (ele_en*ele_en) / ((ele_en-gamma_en)*(ele_en-gamma_en))) * kv(2./3., v)
  term3 = - ele_en / (ele_en-gamma_en) * integrate.quad(K1By3, v, ptive_inf)[0]
  result = term1 * (term2 + term3)
  return result

# --------------------------------------------------------- #
# --------------------------------------------------------- #
# --------------------------------------------------------- #

import os
from os import listdir
from os.path import isfile, join

mypath = os.getcwd()
print('From: ' + mypath)

bin_length = 9

#########################
# EXTRACT LCFA ENERGIES #
#########################

lcfa_nrgs = np.load('./LCFA_tot.npz')['nrgs']

lcfa_points = lcfa_nrgs[0]
lcfa_nrgs = np.delete(lcfa_nrgs, 0)
gamma_lcfa = lcfa_nrgs[0]
lcfa_nrgs = np.delete(lcfa_nrgs, 0)

# HISTOGRAMMING ############

# Binning of the photon emitted gammas read from file
hist, bins = np.histogram(lcfa_nrgs, bins =int(gamma_lcfa/bin_length), range=(0, gamma_lcfa), density= False)

# hist, bins = np.histogram(phtn_gammas, bins =int(gamma/4.995), density= False)
center = (bins[:-1] + bins[1:]) / 2

# zeros = np.where(hist==0)[0]
# print(zeros, zeros.shape)
# for index in zeros:
#     print(center[index])

dbin = center[1] - center[0]
print('bin used to sample the emissions: ', dbin)
norm_track = sum(hist*dbin) # area under the integration region

hist_lcfa = hist / norm_track
center_lcfa = center / gamma_lcfa

##########################
# EXTRACT BLCFA ENERGIES #
##########################

blcfa_nrgs = np.load('./BLCFA_tot.npz')['nrgs']

blcfa_points = blcfa_nrgs[0]
blcfa_nrgs = np.delete(blcfa_nrgs, 0)
gamma_blcfa = blcfa_nrgs[0]
blcfa_nrgs = np.delete(blcfa_nrgs, 0)

#logging
emission_ratio = blcfa_points / lcfa_points
print('LCFA gamma:', gamma_lcfa, 'BLCFA gamma:', gamma_blcfa)
print('#LCFA emissions:', lcfa_points, '#BLCFA emissions:', blcfa_points, 'ratio #BLCFA emissions / #LCFA emissions: ', emission_ratio)

# HISTOGRAMMING ############

# Binning of the photon emitted gammas read from file
hist, bins = np.histogram(blcfa_nrgs, bins =int(gamma_blcfa/bin_length), range=(0, gamma_blcfa), density= False)

# hist, bins = np.histogram(phtn_gammas, bins =int(gamma/4.995), density= False)
center = (bins[:-1] + bins[1:]) / 2

# zeros = np.where(hist==0)[0]
# print(zeros, zeros.shape)
# for index in zeros:
#     print(center[index])

dbin = center[1] - center[0]
print('bin used to sample the emissions: ', dbin)
norm_track = sum(hist*dbin) # area under the integration region

hist_blcfa = hist / norm_track * emission_ratio
center_blcfa = center / gamma_blcfa

#########################
# PLOTTING DISTRIBUTION #
#########################
fig, ax = plt.subplots(figsize = (12, 8))

ax.tick_params(axis='both', direction='in')
plt.xlabel(r'$\omega / \gamma_e$')
plt.ylabel(r'$\frac{dP}{d\omega}$')

ax.plot(center_lcfa, hist_lcfa, '--', color= colors[0],  label='LCFA', linewidth=2)
ax.plot(center_blcfa, hist_blcfa, '--', color = colors[3], label='BLCFA', linewidth=2) 

# ax.set_xlim([3.*bin_length, lcfa_limit + 2.*bin_length])
# ax.set_ylim([0, 0.001])

ax.set_yscale('log')
ax.set_xscale('log')

plt.legend()
plt.tight_layout()

plot_name = 'log_scale_comparison_a0=3.png'
plt.savefig(plot_name)


plt.clf()






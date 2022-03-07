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

file_names = [f for f in listdir(mypath) if (isfile(join(mypath, f)) and f.endswith('.txt'))]

print('searching in', file_names)

x = []
number_of_points = 0

for name in file_names:
    with open(name) as f:
        emitting_gamma = float(f.readline())
        lines = f.readlines()
        #skip empty lines##########
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)
        ###########################
        for line in lines:
            number_of_points = number_of_points + 1
            array = line.split()
            #print(name, line)
            x.append(float(array[0]))

print('particle gamma:', emitting_gamma)
print('number of recorded emissions:', number_of_points)

x = np.array(x)
x = np.insert(x, 0, emitting_gamma)
x = np.insert(x, 0, number_of_points)

#CHANGE THIS NAME ACCORDING TO HOW YOU WANT TO
#SAVE THE npz 
np.savez('BLCFA_tot.npz', nrgs = x)
#np.savez('LCFA_tot.npz', nrgs = x)

x = np.delete(x, 0)
x = np.delete(x, 0)

plt.rc('text', usetex=True)

#############################
# HISTOGRAMMING
#############################

bin_length = 0.9

# Binning of the photon emitted gammas read from file
hist, bins = np.histogram(x, bins =int(emitting_gamma/bin_length), range=(2, emitting_gamma-1), density= False)

# hist, bins = np.histogram(phtn_gammas, bins =int(gamma/4.995), density= False)
center = (bins[:-1] + bins[1:]) / 2

# zeros = np.where(hist==0)[0]
# print(zeros, zeros.shape)
# for index in zeros:
#     print(center[index])

dbin = center[1] - center[0]
norm_track = sum(hist*dbin) # area under the integration region
hist = hist / norm_track
print('bin used to sample the emissions: ', dbin)

#**** Analytical estimate ****
"""

chi = 1.3340

wr_dGamma = []
ph_en_range = center

for ph_en in ph_en_range:
	wr_dGamma.append(wr_dGamma_coeff* Wr_dGamma(ph_en, chi, emitting_gamma)) # gamma is the normalised electron energy

#is this time useful???
# time = 17

wr_dGamma = np.array(wr_dGamma) # s^-1
# dt = time * namelist_dt / namelist_wr # s
# wr_dGamma = wr_dGamma * dt 
dph_en = ph_en_range[1] - ph_en_range[0]
norm = sum(wr_dGamma * dph_en) # area under the integration region
wr_dGamma  = wr_dGamma / norm
"""
#**** ******************* ****

####################
#spectrum
####################
fig, ax = plt.subplots(figsize = (12, 8))

ax.tick_params(axis='both', direction='in')
plt.xlabel(r'$E_{\gamma_s}$ (mc$^2$)')
plt.ylabel(r'$E_{bin} \frac{dN}{dE}$')

ax.plot(center, hist*center, '-', color= colors[0],  label='test result', linewidth=2)
#**************
#ax.plot(ph_en_range, wr_dGamma*ph_en_range, '-', color = colors[2], label='Analytical', linewidth=2) 
#**************
plt.legend()
plt.tight_layout()

plot_name = 'bin=' + str(bin_length) + '_' + str(number_of_points/1000000) + 'M_spec.png'
plt.savefig(plot_name)

####################
#distribution
####################
fig, ax = plt.subplots(figsize = (12, 8))

ax.tick_params(axis='both', direction='in')
plt.xlabel(r'$E_{\gamma_s}$ (mc$^2$)')
plt.ylabel(r'$\frac{dN}{dE}$')

ax.plot(center, hist, '-', color= colors[0],  label='test result', linewidth=2)
#**************
#ax.plot(ph_en_range, wr_dGamma, '-', color = colors[2], label='Analytical', linewidth=2) 
#**************
plt.legend()
plt.tight_layout()

plot_name = 'bin=' + str(bin_length) + '_' + str(number_of_points/1000000) + 'M_distr.png'
plt.savefig(plot_name)


plt.clf()






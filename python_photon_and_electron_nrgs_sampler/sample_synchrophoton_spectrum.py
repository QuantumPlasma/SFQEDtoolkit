import struct
import numpy as np
import matplotlib.pyplot as plt

from scipy.special import kv
from scipy.integrate import quad
import scipy.integrate as integrate

import os

path = os.getcwd()
print('From: ' + path)

bytes_per_double = 8
bytes_per_int = 4

root_folder = './'
file_name_1 = 'electrons'
file_suffix = '.txt'


################## single file reading
x = []
#y = []

number_of_points = 0
with open(root_folder + file_name_1 + file_suffix) as f:
    chi = float(f.readline())
    gamma = float(f.readline())      
    #namelist_dt = float(f.readline())   
    #namelist_wr = float(f.readline())
    lines = f.readlines()
    for line in lines:
    	tmp = float(line)
    	if tmp >= 0 and tmp <= 1:
    	#if True:
            number_of_points = number_of_points + 1
            x.append(tmp)
        #y.append(float(array[1]))
        #z.append(float(array[2]))

e_gammas = np.array(x)

print('number of recorded emissions: ', number_of_points)

#print('x: ')
#for coeff in x_diff:
#    print(coeff)


plt.rc('text', usetex=True)

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
  
  
def k1_3(s):
  return kv(1./3., s)
  
def k2_3(s):
    return kv(2./3., s)
    
def pair_spec(chi, e_nrg):  
    y = 1. / (chi * e_nrg * (1 - e_nrg))
    term1 = integrate.quad(k1_3, (2./3.)*y, ptive_inf)[0]
    term2 = (e_nrg**2 + (1 - e_nrg)**2) / (e_nrg * (1 - e_nrg)) * k2_3( (2./3.)*y )
    return (term1 + term2)

# --------------------------------------------------------- #
# --------------------------------------------------------- #
# --------------------------------------------------------- #

bin_length = 7

# Binning of the photon emitted gammas read from file
hist, bins = np.histogram(e_gammas, bins = 200, range=(0, 1), density= False)

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

# Analytical estimate
wr_dGamma = []
ph_en_range = center

for ph_en in ph_en_range:
	wr_dGamma.append(wr_dGamma_coeff * Wr_dGamma(ph_en, chi, emitting_gamma)) # gamma is the normalised electron energy


wr_dGamma = np.array(wr_dGamma) # s^-1
dph_en = ph_en_range[1] - ph_en_range[0]
norm = sum(wr_dGamma * dph_en) # area under the integration region
wr_dGamma  = wr_dGamma / norm


#plot extrema
x1 = np.amin(center)
x2 = np.amax(center)
print(x1, x2)
ticks_number = 10
dx = (x2 - x1) / ticks_number

#y1 = np.amin(hist)
#y2 = np.amax(hist)
y1 = np.amin(wr_dGamma)
y2 = np.amax(wr_dGamma)
print(y1, y2)
ticks_number = 10
dy = (y2 - y1) / ticks_number


#spectrum
fig, ax = plt.subplots(figsize = (12, 8))

ax.tick_params(axis='both', direction='in')
plt.xlabel(r'$E_{\gamma_s}$ (mc$^2$)')
plt.ylabel(r'$E_{bin} \frac{dN}{dE}$')

ax.plot(center, hist*center, '-', color= colors[0],  label='test result', linewidth=2)
ax.plot(ph_en_range, wr_dGamma*ph_en_range, '-', color = colors[2], label='Analytical', linewidth=2)
plt.xticks(np.arange(x1 - dx, x2 + 2*dx, dx))
plt.yticks(np.arange(y1 - dy, y2 + 2*dy, dy))
plt.legend()
plt.tight_layout()

plot_name = 'chi=' + str(chi) + '_' + str(number_of_points) + '_spec.png'
plt.savefig(plot_name)


#distribution
fig, ax = plt.subplots(figsize = (12, 8))

ax.tick_params(axis='both', direction='in')
plt.xlabel(r'$E_{\gamma_s}$ (mc$^2$)')
plt.ylabel(r'$\frac{dN}{dE}$')

ax.plot(center, hist, '-', color= colors[0],  label='test result', linewidth=2)
ax.plot(ph_en_range, wr_dGamma, '-', color = colors[2], label='Analytical', linewidth=2)
plt.xticks(np.arange(x1 - dx, x2 + 2*dx, dx))
plt.yticks(np.arange(y1 - dy, y2 + 2*dy, dy))
plt.legend()
plt.tight_layout()

plot_name = 'chi=' + str(chi) + '_' + str(number_of_points) + '_distr.png'
plt.savefig(plot_name)

plt.show()

plt.clf()






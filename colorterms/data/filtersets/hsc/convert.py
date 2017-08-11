#!/usr/bin/env python

import numpy as np
import pylab as plt
from scipy import interpolate

""" Utility script to apply various correction to the filter transmission
function :
  - Mirror reflectivity
  - CCD Quantunm efficiency
  - Optical corrector transmission
"""

filters = ['HSC-g', 'HSC-r', 'HSC-i', 'HSC-z', 'HSC-Y']
colors = ['green', 'red', 'firebrick', 'magenta', 'purple']

QE_file = "qe_ccd_HSC.txt"
Mirror_file = "M1-2010s.txt"
Corrector_file = "HSCWFCTx.dat"

# Load auxilliary data
dQE = np.loadtxt(QE_file, unpack=True)
dM = np.loadtxt(Mirror_file, unpack=True)
dC = np.loadtxt(Corrector_file, unpack=True)

# convert to the right set of units (WL in nm and transmission normalized to 1)
dQE[0] = dQE[0]/10.
dC[1] = dC[1]/100.

# plot
fig, (ax0) = plt.subplots(ncols=1, figsize=(5, 5))
ax0.set_xlim([320., 1100])
ax0.set_ylim([0., 1.])
ax0.plot(dQE[0], dQE[1], label = "QE")
ax0.plot(dM[0], dM[1], label = "Mirror")
ax0.plot(dC[0], dC[1], label = "Corrector")
ax0.set_xlabel("Wavelength (nm)")
ax0.set_ylabel("Transmission")
ax0.legend(loc='best')

plt.show()

# Interpolate transmission values
fQE = interpolate.interp1d(dQE[0], dQE[1], fill_value='extrapolate')
fM = interpolate.interp1d(dM[0], dM[1], fill_value='extrapolate')
fC = interpolate.interp1d(dC[0], dC[1], fill_value='extrapolate')

#Aplly to all filters
fig, ax = plt.subplots(ncols=1, figsize=(10, 8))
ax.set_xlim([320., 1200.])
ax.set_ylim([0., 100.])
for i, filt in enumerate(filters):
    fn = filt + "_raw.dat"
    dFilter = np.loadtxt(fn, unpack=True)
    dFilter_new = np.copy(dFilter)
    dFilter_new[1] = dFilter[1]*fQE(dFilter[0])*fM(dFilter[0])*fC(dFilter[0])

    #clean up negative values
    dFilter_new[1][dFilter_new[1] < 0] = 0

    # plot
    ax.plot(dFilter[0], dFilter[1], label = "%s"%filt, c = colors[i])
    ax.plot(dFilter_new[0], dFilter_new[1], c = colors[i])

    # write out the new transmission data filters
    fn = filt + ".dat"
    np.savetxt(fn, np.c_[dFilter_new[0], dFilter_new[1]],  fmt='%.2f  ', delimiter=' ', newline='\n')

ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Transmission (%)")
ax.legend(loc='best')
plt.show()

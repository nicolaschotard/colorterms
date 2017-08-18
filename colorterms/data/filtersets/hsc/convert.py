#!/usr/bin/env python


import numpy as np
import pylab as plt
from scipy import interpolate


"""
Utility script to apply various correction to the filter transmission function.

Corrections:
  - Mirror reflectivity
  - CCD Quantunm efficiency
  - Optical corrector transmission
"""


if __name__ == '__main__':

    filters = ['HSC-g', 'HSC-r', 'HSC-i', 'HSC-z', 'HSC-Y']
    colors = ['green', 'red', 'firebrick', 'magenta', 'purple']

    # Load auxilliary data
    dQE = np.loadtxt("qe_ccd_HSC.txt", unpack=True)  # CCD Quantunm efficiency
    dM = np.loadtxt("M1-2010s.txt", unpack=True)     # Mirror reflectivity
    dC = np.loadtxt("HSCWFCTx.dat", unpack=True)     # Optical corrector transmission

    # convert to the right set of units (WL in nm and transmission normalized to 1)
    dQE[0] = dQE[0]/10.
    dC[1] = dC[1]/100.

    # plot
    fig0, ax0 = plt.subplots(ncols=1, figsize=(9,  7))
    ax0.set_xlim([320., 1100])
    ax0.set_ylim([0., 1.])
    ax0.plot(dQE[0], dQE[1], label = "CCD Quantunm efficiency")
    ax0.plot(dM[0], dM[1], label = "Mirror reflectivity")
    ax0.plot(dC[0], dC[1], label = "Optical corrector transmission")
    ax0.set_xlabel("Wavelength (nm)")
    ax0.set_ylabel("Transmission")
    ax0.legend(loc='best')

    # Interpolate transmission values
    fQE = interpolate.interp1d(dQE[0], dQE[1], fill_value='extrapolate')
    fM = interpolate.interp1d(dM[0], dM[1], fill_value='extrapolate')
    fC = interpolate.interp1d(dC[0], dC[1], fill_value='extrapolate')

    #Aplly to all filters
    fig1, ax1 = plt.subplots(ncols=1, figsize=(9, 7))
    ax1.set_xlim([320., 1200.])
    ax1.set_ylim([0., 100.])
    for i, filt in enumerate(filters):
        fn = filt + "_raw.dat"
        dFilter = np.loadtxt(fn, unpack=True)
        dFilter_new = np.copy(dFilter)
        dFilter_new[1] = dFilter[1] * fQE(dFilter[0]) * fM(dFilter[0]) * fC(dFilter[0])

        #clean up negative values
        dFilter_new[1][dFilter_new[1] < 0] = 0

        # plot
        if i == 0:
            a1, = ax1.plot(dFilter[0], dFilter[1],
                           label="Uncorrected", c='k', lw=0.8)
            a2, = ax1.plot(dFilter_new[0], dFilter_new[1],
                           label='Corrected', ls='--', c='k', lw=0.8)
        ax1.plot(dFilter[0], dFilter[1], c=colors[i])
        ax1.plot(dFilter_new[0], dFilter_new[1], ls='--', c=colors[i])
        ax1.annotate("%s" % filt, (np.mean(dFilter_new[0]), 80), horizontalalignment='center')

        # write out the new transmission data filters
        fn = filt + ".dat"
        np.savetxt(fn, np.c_[dFilter_new[0], dFilter_new[1]],
                   fmt='%.2f  ', delimiter=' ', newline='\n')

        ax1.set_xlabel("Wavelength (nm)")
        ax1.set_ylabel("Transmission (%)")
        ax1.legend(loc='best')
    fig1.savefig("hsc_corrections.png")
    fig1.savefig("hsc_filters.png")
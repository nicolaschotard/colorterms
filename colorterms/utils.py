#!/usr/bin/env python


import numpy as np


class OneSpec(object):

    """
    Single spectrum object with basic atributes.

    lbda: Wavelength
    flux: Flux or transmission
    var:  Variance associated t o the flux (optionnal)
    otype: Type of the object loaded (optionnal)
    """

    def __init__(self, lbda, flux, var=None, otype=None):
        """
        Create a OneSpec instance with wavelength (lbda), flux (flux) and variance (var).
        """
        self.lbda = np.array(lbda, dtype='float')
        self.flux = np.array(flux, dtype='float')
        self.var = var
        self.otype = otype

    def mean_wlength(self):
        """Return effective mean wavelength."""
        return np.sum(self.lbda * self.flux) / np.sum(self.flux)

    def min_wlength(self):
        """Return lower wavelength boundary."""
        return np.min(self.lbda)

    def max_wlength(self):
        """Return higher wavelength boundary."""
        return np.max(self.lbda)


def integ_photons(lbda, flux, step, flbda, filter):

#    if flbda[0] < lbda[0] or flbda[-1] > lbda[-2]:
#        print('Error: %f<%f or %f>%f'%\
#              (flbda[0], lbda[0], flbda[-1], lbda[-2]))
#        return None
    filter_interp = np.interp(lbda, flbda, filter)
    dphotons = (filter_interp * flux) * lbda * 5.006909561e7
    if step is None:
        return np.trapz(dphotons, lbda)
    else:
        return np.sum(dphotons*step)


def integ_photons_variance(lbda, var, step, flbda, filt):

    if flbda[0] < lbda[0] or flbda[-1] > lbda[-2]:
        return None
    filter_interp = np.interp(lbda, flbda, filt)
    dphotons = ((filter_interp * lbda * 5.006909561e7)**2) * var
    return np.sum(dphotons * (step**2))

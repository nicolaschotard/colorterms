#!/usr/bin/env python


import numpy as np


class Spectrum(object):

    """
    Single spectrum object with basic atributes.

    lbda: Wavelength
    flux: Flux or transmission
    var:  Variance associated t o the flux (optionnal)
    object_name: Name of the object (optional)
    object_type: Type of the object (optionnal)
    """

    def __init__(self, lbda, flux, var=None, object_name="", object_type=""):
        """
        Create a Spectrum instance with wavelength (lbda), flux (flux) and variance (var).
        """
        self.lbda = np.array(lbda, dtype='float')
        self.flux = np.array(flux, dtype='float')
        self.var = None if var is None else np.array(var, dtype='float')
        self.steps = self.lbda[1:] - self.lbda[:-1]
        self.constant_step = len(set(self.steps)) == 1
        self.step = None if not self.constant_step else self.steps[0]
        self.object_name = object_name
        self.object_type = object_type

    def mean_wlength(self):
        """Return effective mean wavelength."""
        return np.sum(self.lbda * self.flux) / np.sum(self.flux)

    def min_wlength(self):
        """Return lower wavelength boundary."""
        return np.min(self.lbda)

    def max_wlength(self):
        """Return higher wavelength boundary."""
        return np.max(self.lbda)


class Magnitude(object):
    def __init__(self, spectrum, filters):
        self.spectrum = spectrum
        self.filters = filters
        
        
    def mag(self, step=None, syst='megacam', filt='g'):
        """
        Computes the magnitude of a given spectrum.

        :param array lbda: array of wavelength
        :param array flux: array of flux. Same length as lbda
        :param array var: array of variance. Same length as lbda
        :param array/float int: binning
        :param string syst: the system of filter. See self.filters
        :param string/OneSpec filter: the filter. Must be in the filter set
                                      system, or any OneSpec object.

        :return: magnitude (error)
        """

        # check filter
        filto = self.filters._check_filter(syst, filt)

        photons = integ_photons(self.spectrum.lbda, self.spectrum.flux, step, filto.lbda,
                                filto.flux)
        # refphotons = utils.integ_photons(self.RefSpec.lbda, self.RefSpec.flux,
        #                                  None, filt.lbda, filt.flux)

        if photons is None:  # or refphotons is None:
            if self.spectrum.var is None:
                return -float(np.inf), float(np.inf)
            else:
                return -float(np.inf), float(np.inf)

        outmag = -2.5 / np.log(10) * np.log(photons)  # / refphotons)
        # if self.magoffsets is not None and self.magoffsets.has_key(syst):
        #    if self.magoffsets[syst].has_key(filter):
        #        outmag += self.magoffsets[syst][filter]

        if self.spectrum.var is not None:
            var = integ_photons_variance(self.spectrum.lbda, self.spectrum.var,
                                         self.spetrum.step, filto.lbda, filto.flux)
            magerr = 2.5 / np.log(10) * np.sqrt(var) / photons
            return float(outmag), float(magerr)
        else:
            return float(outmag), None

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

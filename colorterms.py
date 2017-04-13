#!/usr/bin/env python

import os
import numpy as np
import yaml

JOIN = os.path.join

class OneSpec:

    def __init__(self, lbda, flux, var=None):

        self.x = np.array(lbda, dtype='float')
        self.y = np.array(flux, dtype='float')
        self.lbda = self.x
        self.flux = self.y
        self.var = var

    def mean_wlength(self):

        return np.sum(self.lbda * self.flux) / np.sum(self.flux)

class Spectrum(object):

    def __init__(self, lbda, flux, var=None):

        self.lbda = np.array(lbda)
        self.flux = np.array(flux)
        self.var = np.array(var)

    # Filters ============
    def load_filters(self, path):
        """Load all avalaible filter sets"""
        
        # Check first if the filters are already loaded
        if not hasattr(self,'filters'):
            # are the filter set already loaded
            if FILTERS is None:
                # load all avalaible filter sets
                FILTERS = Filters(path)
            self.filters = FILTERS

    def mag(self, var=None, step=None,
            syst='STANDARD', filter='B'):
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
        filt = self._check_filter(syst, filter)

        photons = integ_photons(lbda, flux, step, filt.lbda, filt.flux)
        refphotons = integ_photons(self.RefSpec.lbda, self.RefSpec.flux,
                                   None, filt.lbda, filt.flux)

        if photons is None or refphotons is None:
            if var is None:
                return -float(np.inf), float(np.inf)
            else:
                return -float(np.inf), float(np.inf)

        outmag = -2.5 / np.log(10) * np.log(photons / refphotons)
        if self.magoffsets is not None and self.magoffsets.has_key(syst):
            if self.magoffsets[syst].has_key(filter):
                outmag += self.magoffsets[syst][filter]

        if var is not None:
            var = integ_photons_variance(lbda, var, step, filt.lbda, filt.flux)
            magerr = 2.5 / np.log(10) * np.sqrt(var) / photons
            return float(outmag), float(magerr)
        else:
            return float(outmag), None

    def magerrfit(self, x0, x0err):
        """computes the magnitude error coming from salt2 fit
        (only x0 component is included)"""

        return 2.5 / np.log(10.) * x0err / x0
        
class Filters(object):
        
    def read_filterset_descriptions(path="filtersets"):
        """
        Get all filter sets description
        
        Return a dictionnary containing the different sets of filters
        """
        self.filtersets = yaml.load(path + "/description.yaml")
        
     def load_filters(self, path):
        """
        Load the filter names and transmission stored in the model directory.
        Also load the reference spectrum used in the model (for magnitude).
        See the following attributes:
        - FilterWheels
        - Filters
        - RefSpec
        """
        
        # Check first if the filters are already loaded
        if hasattr(self,'filters'):
            return

        # will load all the avalaible filters
        self.filters = load_filter_sets(self.dir, self.cards)
        
        # Load the reference spectrum:
        lbda, flux = np.loadtxt(JOIN(self.dir, self.cards['REF_SPECTRUM']),
                                unpack=True)
        self.RefSpec = OneSpec
        (lbda,flux)   
        
    def _check_filter(self, syst, filt):
        """
        Make sure that a given filter is a OneSpec object or is in the list
        of already defined filters.
        """

        def check_attributes(f):

            return hasattr(f,'mean_wlength') & \
                   hasattr(f,'lbda') & \
                   hasattr(f,'flux')

        # if filt is a OneSpec object
        if check_attributes(filt): 
            return filt
        
        # check the system of filters
        if not self.filters.has_key(syst):
            raise KeyError("Selected system (%s) not in self.filters"%syst)
        
        # Check if the asked filter is in the list
        if self.filters[syst].has_key(filt):
            return self.filters[syst][filt]
        else:
            raise KeyError("Selected filter (%s) does not exist or has the '\
                           'wrongattributes"%filt)

# Magnitudes utilities====================

def integ_photons(lbda, flux, step, flbda, filter):

    if flbda[0] < lbda[0] or flbda[-1] > lbda[-2]:
        print 'Error: %f<%f or %f>%f'%\
              (flbda[0], lbda[0], flbda[-1], lbda[-2])
        return None
    filter_interp = np.interp(lbda, flbda, filter)
    dphotons = (filter_interp * flux) * lbda * 5.006909561e7
    if step is None:
        return np.trapz(dphotons,lbda)
    else:
        return np.sum(dphotons*step)

def integ_photons_variance(lbda, var, step, flbda, filter):

    if flbda[0] < lbda[0] or flbda[-1] > lbda[-2]:
        return None
    filter_interp = np.interp(lbda, flbda, filter)
    dphotons = ((filter_interp * lbda * 5.006909561e7)**2) * var
    return np.sum(dphotons * (step**2))

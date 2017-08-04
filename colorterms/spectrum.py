#!/usr/bin/env python


import numpy as np
from . import utils


class Spectrum(object):
    """Create a spectrum instance with a few useful methods"""

    filters = None

    def __init__(self, lbda, flux, var=None, fpath='filtersets'):

        self.lbda = np.array(lbda)
        self.flux = np.array(flux)
        self.var = np.array(var)
        self._load_filters(fpath)

    def _load_filters(self, path):
        """Load all avalaible filter sets"""

        # Check first if the filters are already loaded
        if self.filters is None:
            Spectrum.filters = Filters(path_to_filters=path, verbose=False)

    def mag(self, var=None, step=None, syst='megacam', filt='g'):
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

        photons = utils.integ_photons(self.lbda, self.flux, step, filto.lbda,
                                      filto.flux)
        # refphotons = utils.integ_photons(self.RefSpec.lbda, self.RefSpec.flux,
        #                                  None, filt.lbda, filt.flux)

        if photons is None:  # or refphotons is None:
            if var is None:
                return -float(np.inf), float(np.inf)
            else:
                return -float(np.inf), float(np.inf)

        outmag = -2.5 / np.log(10) * np.log(photons)  # / refphotons)
        # if self.magoffsets is not None and self.magoffsets.has_key(syst):
        #    if self.magoffsets[syst].has_key(filter):
        #        outmag += self.magoffsets[syst][filter]

        if var is not None:
            var = utils.integ_photons_variance(lbda, var, step, filt.lbda, filt.flux)
            magerr = 2.5 / np.log(10) * np.sqrt(var) / photons
            return float(outmag), float(magerr)
        else:
            return float(outmag), None
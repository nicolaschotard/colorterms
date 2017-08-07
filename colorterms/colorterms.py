#!/usr/bin/env python


import numpy as np
import pylab as plt
from . import spectools


class Colorterms(object):

    def __init__(self, catalogs, filters):
        """Initialization."""
        # Initiale inputs
        self.catalogs = [catalogs] if not isinstance(catalogs, list) else catalogs
        self.filters = filters

        # Create mags object for all spectra of all catalogs
        self.mag_catalogs = [[spectools.Magnitude(spec, filters) for spec in cat.spectra]
                             for cat in self.catalogs]

        # Compute magnitudes for all spectra of all catalogs, and for all systems and filters
        self.magnitudes = {}
        for i, mag_catalog in enumerate(self.mag_catalogs):
            cmag = self.magnitudes[self.catalogs[i].name] = {}
            for syst in filters.filters:
                cmag[syst] = {}
                for filt in filters.filters[syst]:
                    cmag[syst][filt] = np.array([mag.mag(syst=syst, filt=filt)[0]
                                                 for mag in mag_catalog])

    def _pair_filters_between_systems(self, first_filterset, second_filterset):
        """Pair filters from one system to an other."""
        # mean_wlength should be close enough. Diff < 100A ?
        # range should also be close enough.
        filters = self.filters.filters
        for filt in filters[first_filterset]:
            pass

    def _pair_filters_for_colors(self, first_filterset, second_filterset):
        """Get filter pairs for color in a given system."""
        pass

    def compute_colorterms(first_filterset, second_filterset):
        self.first_filterset = first_filterset
        self.second_filterset = second_filterset

    def plot_c_vs_magdiff(first_filterset, second_filterset):
        pass

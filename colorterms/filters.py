#!/usr/bin/env python


import numpy as np
import yaml
import pylab as plt
from . import utils


class Filters(object):

    def __init__(self, path_to_filters="filtersets", verbose=True):
        """Load all available filter sets."""
        self.path_to_filters = path_to_filters
        self.load_filters(verbose=verbose)
        self.ordered = self.order_by_wlength()

    def _read_filterset_descriptions(self):
        """
        Get all filter sets description.

        Return a dictionnary containing the different sets of filters
        """
        self.filtersets = yaml.load(open(self.path_to_filters + "/description.yaml"))

    def load_filters(self, verbose=True):
        """
        Load the filter names and transmission stored in the model directory.
        Also load the reference spectrum used in the model (for magnitude).
        See the following attributes:
        - FilterWheels
        - Filters
        - RefSpec
        """
        # Get the filter description
        self._read_filterset_descriptions()

        # Check first if the filters are already loaded
        if hasattr(self, 'filters'):
            return

        # will load all the avalaible filters
        self.filters = {fset: self._load_filter_set(fset, verbose=verbose)
                        for fset in self.filtersets}

    def _load_filter_set(self, fset, verbose=True):
        """Load a given filter set."""
        if verbose:
            print("INFO: Loading %s filter set" % fset)
        data = {}
        for filt in self.filtersets[fset]:
            if verbose:
                print(" - loading %s" % filt)
            d = np.loadtxt("%s/%s/%s" % (self.path_to_filters, fset,
                                         self.filtersets[fset][filt]), unpack=True)
            data[filt] = utils.OneSpec(d[0], d[1])
        return data

    def _check_filter(self, syst, filt):
        """
        Make sure that a given filter is a OneSpec object or is in the list
        of already defined filters.
        """

        def check_attributes(f):
            """Check attributes of the given filter set."""
            return hasattr(f, 'mean_wlength') & \
                hasattr(f, 'lbda') & \
                hasattr(f, 'flux')

        # if filt is a OneSpec object
        if check_attributes(filt):
            return filt

        # check the system of filters
        if syst not in self.filters:
            raise KeyError("Selected system (%s) not in self.filters"%syst)

        # Check if the asked filter is in the list
        if filt in self.filters[syst]:
            return self.filters[syst][filt]
        else:
            raise KeyError("Selected filter (%s) does not exist or has the '\
                           'wrongattributes"%filt)

    def plot_filters(self):
        """Simple transmission plots."""
        for syst in self.filters:
            fig = plt.figure(dpi=150)
            ax = fig.add_subplot(111,
                                 xlabel='Wavelenght',
                                 ylabel='Transmission',
                                 title=syst)
            for filt in self.filters[syst]:
                ax.plot(self.filters[syst][filt].lbda,
                        self.filters[syst][filt].flux,
                        label=filt)
                ax.legend(loc='best')

    def order_by_wlength(self):
        """
        Order filters by average wavelength values for each system.

        Return a dictionnary of orderer arrays.
        """
        ofilt = {}
        for syst in self.filters:
            filters = np.array(sorted(self.filters[syst]))
            meanw = np.array([self.filters[syst][filt].mean_wlength() for filt in filters])
            ofilt[syst] = filters[np.argsort(meanw)]
        return ofilt

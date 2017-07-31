#!/usr/bin/env python


from glob import glob
import numpy as np
import yaml


CATALOGS = {}


class OneSpec(object):

    def __init__(self, lbda, flux, var=None):

        self.x = np.array(lbda, dtype='float')
        self.y = np.array(flux, dtype='float')
        self.lbda = self.x
        self.flux = self.y
        self.var = var
        self.object_type = None

    def mean_wlength(self):

        return np.sum(self.lbda * self.flux) / np.sum(self.flux)


class Spectrum(object):

    def __init__(self, lbda, flux, var=None, fpath='filtersets'):

        self.lbda = np.array(lbda)
        self.flux = np.array(flux)
        self.var = np.array(var)
        self._load_filters(fpath)

    def _load_filters(self, path):
        """Load all avalaible filter sets"""

        # Check first if the filters are already loaded
        if not hasattr(self, 'filters'):
            # are the filter set already loaded
            self.filters = Filters(path_to_filters=path, verbose=False)

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

        photons = integ_photons(self.lbda, self.flux, step, filto.lbda,
                                filto.flux)
        # refphotons = integ_photons(self.RefSpec.lbda, self.RefSpec.flux,
        #                           None, filt.lbda, filt.flux)

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

    def __init__(self, path_to_filters="filtersets", verbose=True):
        """Load all available filter sets."""
        self.path_to_filters = path_to_filters
        self._load_filters(verbose=verbose)

    def _read_filterset_descriptions(self):
        """
        Get all filter sets description.

        Return a dictionnary containing the different sets of filters
        """
        self.filtersets = yaml.load(open(self.path_to_filters + "/description.yaml"))

    def _load_filters(self, verbose=True):
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
            data[filt] = OneSpec(d[0], d[1])
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

# Magnitudes utilities====================


def integ_photons(lbda, flux, step, flbda, filter):

    if flbda[0] < lbda[0] or flbda[-1] > lbda[-2]:
        print('Error: %f<%f or %f>%f'%\
              (flbda[0], lbda[0], flbda[-1], lbda[-2]))
        return None
    filter_interp = np.interp(lbda, flbda, filter)
    dphotons = (filter_interp * flux) * lbda * 5.006909561e7
    if step is None:
        return np.trapz(dphotons, lbda)
    else:
        return np.sum(dphotons*step)


def integ_photons_variance(lbda, var, step, flbda, filter):

    if flbda[0] < lbda[0] or flbda[-1] > lbda[-2]:
        return None
    filter_interp = np.interp(lbda, flbda, filter)
    dphotons = ((filter_interp * lbda * 5.006909561e7)**2) * var
    return np.sum(dphotons * (step**2))


def load_catalogs(path="catalogs", catalog="gunnstryker", ext=".ascii",
                  desc="lbd,flux"):

    if catalog == "gunnstryker":
        # paths
        spectra = sorted(glob(path + "/" + catalog + "/gs_*.ascii"))
        CATALOGS[catalog] = {int(sp.split('_')[1].split('.')[0]):
                             {'path': sp, 'type': 'unknown',
                              'lbda': None, 'flux': None, 'spec': None}
                             for sp in spectra}

        # spectrum type is any
        spectype = np.loadtxt(path + "/" + catalog + "/gsspectype.ascii",
                              dtype='str', unpack=True)
        for i, gs in enumerate(spectype[0]):
            CATALOGS[catalog][int(gs.split('_')[1].split('.')[0])]['type'] = spectype[1][i]

        # data
        filtersets = '/'.join(path.split('/')[:-1] + ['filtersets'])
        for sp in spectra:
            x, y = np.loadtxt(sp, dtype='float', unpack=True)
            num = int(sp.split('_')[1].split('.')[0])
            CATALOGS[catalog][num]['lbda'] = x
            CATALOGS[catalog][num]['flux'] = y
            CATALOGS[catalog][num]['spec'] = Spectrum(x, y, fpath=filtersets)
    else:  # use input args
        spectra = glob(path + "/" + catalog + "/*.ext")

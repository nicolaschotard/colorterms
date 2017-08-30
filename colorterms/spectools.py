#!/usr/bin/env python


import numpy as np
import pylab as pl


class Spectrum(object):

    """
    Single spectrum object with basic atributes.

    lbda: Wavelength
    flux: Flux or transmission
    var:  Variance associated t o the flux (optionnal)
    object_name: Name of the object (optional)
    object_type: Type of the object (optionnal)
    """

    def __init__(self, lbda, flux, **kwargs):
        """
        Create a Spectrum instance with wavelength (lbda), flux (flux) and variance (var).

        Available kwargs are:
        - var: Variance associated to the spectrum flux
        - object_name: Name of the objet loaded
        - object_type: Type of the object loaded
        """
        var = kwargs.get('var', None)
        self.lbda = np.array(lbda, dtype='float')
        self.flux = np.array(flux, dtype='float')
        self.var = None if var is None else np.array(var, dtype='float')
        self.steps = self.lbda[1:] - self.lbda[:-1]
        self.constant_step = len(set(self.steps)) == 1
        self.step = None if not self.constant_step else self.steps[0]
        self.object_name = kwargs.get("object_name", None)
        self.object_type = kwargs.get("object_type", None)

    def mean_wlength(self):
        """Return effective mean wavelength."""
        return np.sum(self.lbda * self.flux) / np.sum(self.flux)

    def min_wlength(self):
        """Return lower wavelength boundary."""
        return np.min(self.lbda)

    def max_wlength(self):
        """Return higher wavelength boundary."""
        return np.max(self.lbda)

    def fwhm_wlength(self):
        """
        Return the Full-Width Half-Maximum.
        
        Can be used as a proxy of the width of a bandpass filter.
        """
        wlength_range = self.lbda[self.flux > max(self.flux) / 2]
        return wlength_range[-1] - wlength_range[0]


class Magnitude(object):

    """Compute magnitudes for a given spectrum a filter sets."""

    def __init__(self, spectrum, filtersets):
        """
        Initialize a Magnitude object, ready to compute magnitudes for a 
        given spectrum a filter sets.

        :param object spectrum: A spectools.Spectrum object
        :param object filters: A filtersets.Filters object
        """
        self.spectrum = spectrum
        self.filters = filtersets

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
        filto = self.filters.check_filter(syst, filt)

        photons = integ_photons(self.spectrum.lbda, self.spectrum.flux, step,
                                filto.lbda, filto.flux)
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


class Magnitudes(object):
    
    def __init__(self, catalogs, filters):

        self.catalogs = catalogs
        self.filters = filters
        self._mag_catalogs = {}
        self.magnitudes = {}
    

    def compute_magnitudes(self, filtersets=None, catalogs=None):
        """Compute the magnitudes for a given system and catalog.
        
        If 'filters' is None, magnitudes will be computed for all available filter sets.
        If 'catalog_list' is None, magnitudes will be computed for all available catalogs.
        Magnitudes are stored in self.magnitudes.
        """
        catalogs = list(catalogs) if catalogs is not None else self.catalogs.keys()
        # Create mags object for all spectra of all catalogs
        for cat in catalogs:
            if cat in self._mag_catalogs:
                continue
            else:
                self._mag_catalogs[cat] = [Magnitude(spec, self.filters)
                                           for spec in self.catalogs[cat].spectra]

        # Compute magnitudes for all spectra for the input list of catalogs
        # and for all filters of the two input systems
        filtersets = list(filtersets) if filtersets is not None else self.filters.filters.keys()        
        for cat in catalogs:
            print("INFO: Computing magnitudes for the %s catalog" % cat)
            if cat not in self.magnitudes.keys():
                cmag = self.magnitudes[cat] = {}
            else:
                cmag = self.magnitudes[cat]
            for syst in filtersets:
                if syst in cmag:
                    continue
                print(" -> filter set: %s" % syst)
                cmag[syst] = {}
                for filt in self.filters.filters[syst]:
                    cmag[syst][filt] = np.array([mag.mag(syst=syst, filt=filt)[0]
                                                 for mag in self._mag_catalogs[cat]])

    def hists(self, filtersets=None, catalogs=None):
        catalogs = list(catalogs) if catalogs is not None else self.catalogs.keys()
        filtersets = list(filtersets) if filtersets is not None else self.filters.filters.keys()
        for fset in filtersets:
            for filt in self.filters.filters[fset].keys():
                mags, labels = [], []
                for cat in catalogs:
                    if cat not in self.magnitudes.keys():
                        continue
                    elif fset not in self.magnitudes[cat].keys():
                        continue
                    else:
                        mags.append(self.magnitudes[cat][fset][filt])
                        labels.append(cat)
                if len(mags) != 0:
                    fig = pl.figure()
                    ax = fig.add_subplot(111)
                    ax.xlabel('%s - %s magnitudes' % (fset, filt))
                    for mag, label in zip(mags, labels):
                        ax.hist(mag, label=label)
                    ax.legend(loc='best')
        pl.show()


def integ_photons(lbda, flux, step, flbda, filt):

    # if flbda[0] < lbda[0] or flbda[-1] > lbda[-2]:
    #     print('Error: %f<%f or %f>%f'%\
    #           (flbda[0], lbda[0], flbda[-1], lbda[-2]))
    #     return None
    filter_interp = np.interp(lbda, flbda, filt)
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

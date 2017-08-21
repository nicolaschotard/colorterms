#!/usr/bin/env python


import numpy as np
import pylab as plt
from scipy import polyfit, polyval
from . import spectools


class Colorterms(object):

    def __init__(self, catalogs, filters):
        """Initialization."""
        # Initiale inputs
        self.catalogs = catalogs
        self.filters = filters

        # Initialization of variables
        self.mag_catalogs = {}
        self.magnitudes = {}
        self.colorterms = {}
        self.pairs = {}

    def _compute_magnitudes(self, first_fset, second_fset, catalog_list=None):
        """Compute the magnitudes for the two given system.

        If the 'catalog_list' is not given, magnitudes will be computed for all available catalogs.
        Magnitudes are stored in self.magnitudes.
        """
        # Create mags object for all spectra of all catalogs
        catalogs = catalog_list if catalog_list is not None else self.catalogs.keys()
        for cat in catalogs:
            if cat in self.mag_catalogs:
                continue
            else:
                self.mag_catalogs[cat] = [spectools.Magnitude(spec, self.filters)
                                          for spec in self.catalogs[cat].spectra]

        # Compute magnitudes for all spectra for the input list of catalogs
        # and for all filters of the two input systems
        for cat in catalogs:
            if cat not in self.magnitudes.keys():
                cmag = self.magnitudes[cat] = {}
            else:
                cmag = self.magnitudes[cat]
            for syst in [first_fset, second_fset]:
                if syst in cmag:
                    continue
                cmag[syst] = {}
                for filt in self.filters.filters[syst]:
                    cmag[syst][filt] = np.array([mag.mag(syst=syst, filt=filt)[0]
                                                 for mag in self.mag_catalogs[cat]])

    def _make_pairing(self, first_fset, second_fset):
        """
        Pair filters from one system to an other.

        The mean wavelength of two filters should be close enough. Diff < 100A ?
        Their range should also of the same order.
        """
        if second_fset not in self.pairs:
            self.pairs[second_fset] = {}
        if first_fset not in self.pairs[second_fset]:
            results = {}
            for filt_2 in self.filters.ordered[second_fset]:
                # Mean wavelength for each filter of the first set
                means = np.array([self.filters.filters[first_fset][filt_1].mean_wlength()
                                  for filt_1 in self.filters.ordered[first_fset]])
                # Wavelength difference of the current filter to all filters of the first set
                diff = np.abs(self.filters.filters[first_fset][filt_2].mean_wlength() - means)
                # Argument of the closest one
                amin = np.argmin(diff)
                # Closest one
                closest_filt = self.filters.ordered[first_fset][amin]
                results[filt_2] = {'filter': closest_filt}
                # Colors corresponding to this closest filter
                filts = list(self.filters.ordered[first_fset])
                colors = np.concatenate([((filt, (filts + filts[:1])[i - 1]),
                                          (filt, (filts + filts[:1])[i + 1]))
                                         for i, filt in enumerate(filts)])[1:-1]
                results[filt_2]['colors'] = colors
            self.pairs[second_fset][first_fset] = results

    def _get_data(self, first_fset, second_fset, filt, color, catalogs, cuts):
        """Return valid data."""
        # Short cut
        localdic = self.pairs[second_fset][first_fset][filt]

        # Get the data for all catalogs
        m0 = np.concatenate([self.magnitudes[catalog][second_fset][filt]
                             for catalog in catalogs])
        m1 = np.concatenate([self.magnitudes[catalog][first_fset][localdic['filter']]
                             for catalog in catalogs])
        col = np.concatenate([self.magnitudes[catalog][first_fset][color[0]] - \
                                      self.magnitudes[catalog][first_fset][color[1]]
                              for catalog in catalogs])

        # first make sure to remove all possible inf or None values
        mask = np.isfinite(m0) & np.isfinite(m1) & np.isfinite(col)
        m0, m1, col = m0[mask], m1[mask], col[mask]

        # Apply filters if any
        if cuts is not None:
            mask = self._get_mask(second_fset, filt, m0, cuts)
            mask &= self._get_mask(first_fset, localdic['filter'], m1, cuts)
            mask &= self._get_mask(first_fset, "%s-%s" % (color[0], color[1]), col, cuts)
            m0, m1, col = m0[mask], m1[mask], col[mask]

        return m0, m1, col

    def _get_mask(self, fset, param, data, cuts):
        """
        Apply cuts on a data set.

        fset: The filter set
        param: A filter name of a color name (e.g., g-i)
        data: a numpy array of data to filter
        cuts: a dictionnary of cuts

        return a mask of the same dimension as data
        """
        mask = np.ones(len(data), dtype='bool')
        if fset not in cuts:
            return mask
        if param not in cuts[fset]:
            return mask
        if 'min' in cuts[fset][param]:
            mask &= data >= cuts[fset][param]['min']
        if 'max' in cuts[fset][param]:
            mask &= data <= cuts[fset][param]['max']
        return mask

    def compute_colorterms(self, first_fset, second_fset, catalogs=None, cuts=None):
        """Compute colorterm slopes to go from the first filterset to the second one.

        Transformations are of the following types:
        s1(f) = s2(f') + a * [s2(f') - s2(f'')]

        Where s1(f) and s2(f') must be as close as possible from each other in term of mean
        wavelength and wavelength range. s2(f') and s2(f'') will be side by side filters.

        e.g.: SDSS(g) = MEGACAM(g) + a * [MEGACAM(g) - MEGACAM(r)]

        catalogs: List of catalog to use in the fits
        cuts: A dictionnary containing a possible list of cuts for each filter or color.
              This dictionnary is of the following form, with no mandaotory keys:
              cuts = {'megacam': {'g': {'min': 10, 'max': 22}       # for an individual filter
                                  'g-r': {'min': 0.1, 'max': 1.5}   # for a color
                                  }
                      }
              Filter sets and filters must exist. Colors are of the form 'f1-f2'.
        """
        catalogs = self.catalogs.keys() if catalogs is None else catalogs
        self._make_pairing(first_fset, second_fset)
        self._compute_magnitudes(first_fset, second_fset)
        print("INFO: Computing colorterms to go from %s to %s" % (first_fset, second_fset))
        for filt in self.pairs[second_fset][first_fset]:
            localdic = self.pairs[second_fset][first_fset][filt]
            for color in localdic['colors']:
                print(" FITTING: %s(%s) - %s(%s) = f(%s(%s) - %s(%s))" %
                      (second_fset, filt, first_fset, localdic['filter'],
                       first_fset, color[0], first_fset, color[1]))
                m0, m1, col = self._get_data(first_fset, second_fset, filt, color, catalogs, cuts)
                colfit = Colorfit(m0 - m1, col,
                                  xlabel="%s(%s) - %s(%s)" % (first_fset, color[0],
                                                              first_fset, color[1]),
                                  ylabel="%s(%s) - %s(%s)" % (second_fset, filt, first_fset,
                                                              localdic['filter']),
                                  title="%s, %s filter" % (second_fset, filt))
                colfit.polyfits()
                colfit.plots()
                for order in colfit.polyfits_outputs:
                    print("Order =", order, colfit.polyfits_outputs[order]['params'])


    def plot_magdiff_vs_c(self, first_fset, second_fset, catalogs=None):
        """Magnitude difference as a function of color. DEPRECATED FOR NOW"""
        catalogs = self.catalogs.keys() if catalogs is None else catalogs

        self._make_pairing(first_fset, second_fset)
        self._compute_magnitudes(first_fset, second_fset)
        for filt in self.pairs[second_fset][first_fset]:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            localdic = self.pairs[second_fset][first_fset][filt]
            ax.set_xlabel("%s(%s) - %s(%s)" %
                          (first_fset, localdic['filter'], first_fset, localdic['colors'][0]))
            ax.set_ylabel("%s(%s) - %s(%s)" %
                          (second_fset, filt, first_fset, localdic['filter']))
            ax.set_title("%s, %s filter" %
                         (second_fset, filt))
            for catalog in catalogs:
                x = self.magnitudes[catalog][first_fset][localdic['filter']] - \
                    self.magnitudes[catalog][first_fset][localdic['colors'][0]]
                y = self.magnitudes[catalog][second_fset][filt] - \
                    self.magnitudes[catalog][first_fset][localdic['filter']]
                ax.plot(x, y, 'o', label='%s (%i)' % (catalog, len(x)))
            ax.legend(loc='best')
        plt.show()


class Colorfit(object):

    def __init__(self, magdiff, color, **kwargs):
        """
        Polynomial fit of magnitude difference versus color (i.e, x vs y)

        Possible kwargs:
        xlabel: Label of the `color` argument for plot purpose
        ylabel: Label of the `magdiff` argument for plot purpose
        title: A title for the figure
        """
        self.magdiff = magdiff
        self.color = color
        self.kwargs = kwargs
        self.params = {}
        self.polyfits_outputs = {}

    def polyfits(self, orders="1,2,3"):
        """
        Simple polynomial fits of order 1, 2, 3.

        orders: A integer, a list of string or integer, or a string.
                e.g.: orders = 1 or "1" or "1,2,3" or ["1", "2"] or [1, 2]

        Results saved in self.polyfits_outputs
        """
        if isinstance(orders, str):
            orders = [int(order) for order in orders.split(",")]
        elif isinstance(orders, list):
            orders = [int(order) for order in orders]
        elif isinstance(orders, int):
            pass
        else:
            raise IOError("The 'orders' argument must be a integer, a string or a list")
        for order in orders:
            output = self.polyfits_outputs[order] = {}
            output['params'] = polyfit(self.color, self.magdiff, order)
            output['ymodel'] = polyval(output["params"], self.color)
            output['yresiduals'] = self.magdiff - output['ymodel']
            output['yresiduals_mean'] = np.mean(output['yresiduals'])
            output['yresiduals_std'] = np.std(output['yresiduals'])

    def plots(self):
        """Plot the polynomial fit results."""
        if len(self.polyfits_outputs) == 0:
            raise "INFO: You must run the polyfits method first"
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(self.kwargs.get("xlabel", ""))
        ax.set_ylabel(self.kwargs.get("ylabel", ""))
        ax.set_title(self.kwargs.get("title", "") +
                     ", %i data points" % len(self.color))
        ax.plot(self.color, self.magdiff, 'ok')
        for order in self.polyfits_outputs:
            xsorted = np.linspace(min(self.color), max(self.color), 200)
            ysorted = polyval(self.polyfits_outputs[order]["params"], xsorted)
            ax.plot(xsorted, ysorted, label="order=%i, std=%.3f" %
                    (order, self.polyfits_outputs[order]["yresiduals_std"]))
        ax.legend(loc='best')
        plt.show()

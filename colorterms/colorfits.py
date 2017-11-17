#!/usr/bin/env python


import os
import copy
import pickle
import yaml
import numpy as np
import pylab as plt
from scipy import polyfit, polyval
from . import spectools


class Colorterms(object):

    """Compute color terms to go from one filter set to an other."""

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
        self.results = {}

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
            paired = {}
            for filt_2 in self.filters.ordered[second_fset]:
                # Mean wavelength for each filter of the first set
                means = np.array([self.filters.filters[first_fset][filt_1].mean_wlength()
                                  for filt_1 in self.filters.ordered[first_fset]])
                # Wavelength difference of the current filter to all filters of the first set
                diff = np.abs(self.filters.filters[second_fset][filt_2].mean_wlength() - means)
                # Argument of the closest one
                amin = np.argmin(diff)
                # Closest one
                closest_filt = self.filters.ordered[first_fset][amin]
                paired[str(filt_2)] = {'filter': str(closest_filt)}
                # Colors corresponding to this closest filter
                filts = list(self.filters.ordered[first_fset])
                colors = np.concatenate([((filt, (filts + filts[:1])[i - 1]),
                                          (filt, (filts + filts[:1])[i + 1]))
                                         for i, filt in enumerate(filts)])[1:-1]
                # Only keep colors where the closest filter is also part of the color definition
                colors = [list(c) for c in colors if c[0] == closest_filt]
                paired[str(filt_2)]['colors'] = colors
                self.pairs[second_fset][first_fset] = paired

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
            mask = self._get_mask("%s(%s)" % (second_fset, filt), m0, cuts)
            mask &= self._get_mask("%s(%s)" % (first_fset, localdic['filter']), m1, cuts)
            mask &= self._get_mask("%s(%s) - %s(%s)" % (second_fset, filt,
                                                        first_fset, localdic['filter']),
                                   m0 - m1, cuts)
            mask &= self._get_mask("%s(%s) - %s(%s)" % (first_fset, color[0], first_fset, color[1]),
                                   col, cuts)
            m0, m1, col = m0[mask], m1[mask], col[mask]

        return m0, m1, col

    def _get_mask(self, param, data, cuts):
        """
        Apply cuts on a data set.

        fset: The filter set
        param: A filter name of a color name (e.g., g-i)
        data: a numpy array of data to filter
        cuts: a dictionnary of cuts

        return a mask of the same dimension as data
        """
        mask = np.ones(len(data), dtype='bool')
        iparam = ' - '.join(param.split(' - ')[::-1])  # (a - b) -> (b - a)
        if param in cuts:
            if 'min' in cuts[param]:
                mask &= data >= cuts[param]['min']
            if 'max' in cuts[param]:
                mask &= data <= cuts[param]['max']
        elif iparam in cuts:
            if 'min' in cuts[iparam]:
                mask &= data <= -cuts[iparam]['min']
            if 'max' in cuts[iparam]:
                mask &= data >= -cuts[iparam]['max']
        return mask

    def compute_colorterms(self, first_fset, second_fset, catalogs=None,
                           cuts=None, sigma_clip=None, verbose=False):
        """Compute colorterm slopes to go from the first filterset to the second one.

        Transformations are of the following types:
        s1(f) = s2(f') + a * [s2(f') - s2(f'')]

        Where s1(f) and s2(f') must be as close as possible from each other in term of mean
        wavelength and wavelength range. s2(f') and s2(f'') will be side by side filters.

        e.g.: SDSS(g) = MEGACAM(g) + a * [MEGACAM(g) - MEGACAM(r)]

        catalogs: List of catalog to use in the fits
        cuts: A dictionnary containing a possible list of cuts for each filter or color.
              This dictionnary is of the following form, with no mandatory key:
              cuts = {'megacam(g) - sdss(g)': {'min': 10, 'max': 22}
                      'sdss(g) - sdss(r)': {'min': 0.1, 'max': 1.5}
                     }
              All cuta are of the form: system(filt) - system(filt).
              If 'sdss(g) - megacam(r)' is defined, 'megacam(r) - sdss(g)' will automatically works.
        sigma_clip: Sigma clipping value while fitting for color terms
        """
        catalogs = list(self.catalogs.keys()) if catalogs is None else catalogs
        self._make_pairing(first_fset, second_fset)
        self._compute_magnitudes(first_fset, second_fset)
        print("INFO: Computing colorterms to go from %s to %s" % (first_fset, second_fset))
        for filt in self.pairs[second_fset][first_fset]:
            localdic = self.pairs[second_fset][first_fset][filt]
            localdic['results'] = {}
            for color in localdic['colors']:
                if verbose:
                    print(" FITTING: %s(%s) - %s(%s) = f(%s(%s) - %s(%s))" %
                          (second_fset, filt, first_fset, localdic['filter'],
                           first_fset, color[0], first_fset, color[1]))
                m0, m1, col = self._get_data(first_fset, second_fset, filt, color, catalogs, cuts)
                if not any(m0-m1):
                    print("WARNING: Skipping  these two identical filters: %s(%s) & %s(%s)" % \
                          (second_fset, filt, first_fset, localdic['filter']))
                    continue
                colfit = Colorfit(m0 - m1, col,
                                  xlabel="%s(%s) - %s(%s)" % (first_fset, color[0],
                                                              first_fset, color[1]),
                                  ylabel="%s(%s) - %s(%s)" % (second_fset, filt, first_fset,
                                                              localdic['filter']),
                                  title="%s, %s filter" % (second_fset, filt))
                colfit.polyfits(sigma_clip=sigma_clip)
                # By Catalog Data (for the final plots)
                bcd = np.transpose([self._get_data(first_fset, second_fset, filt, color, [catalog],
                                                   cuts) for catalog in catalogs])
                colfit.plots(dirname="%s_%s" % (first_fset, second_fset),
                             bycat_data=np.transpose([bcd[0], bcd[1], bcd[2], catalogs]))
                results = localdic['results'][",".join(color)] = {}
                for order in colfit.polyfits_outputs:
                    if verbose:
                        print("Order =", order, colfit.polyfits_outputs[order]['params'],
                              " (STD=%.3f)" % colfit.polyfits_outputs[order]['yresiduals_std'])
                    results[order] = colfit.polyfits_outputs[order]
        

        # Order results by rms and print them
        self._order_by_rms(first_fset, second_fset)

    def _order_by_rms(self, first_fset, second_fset):
        
        # Order them by best RMS for each pair
        for filt in self.pairs[second_fset][first_fset]:
            localdic = self.pairs[second_fset][first_fset][filt]
            print(" BEST FIT FOR: %s(%s) - %s(%s) = f(%s(??) - %s(??))" %
                  (second_fset, filt, first_fset, localdic['filter'],
                   first_fset, first_fset))
            colors, results = [], []
            for color in localdic['results']:

                for order in localdic['results'][color]:
                    colors.append("f(%s(%s) - %s(%s)) [%i]" % (first_fset, color.split(',')[0],
                                                               first_fset, color.split(',')[1],
                                                               int(order)))
                    results.append(localdic['results'][color][order]['yresiduals_std'])
            for c, r in zip(np.array(colors)[np.argsort(results)], np.sort(results)):
                print(c, ": RMS=%.3f" % r)

    def build_colorterms_dict(self):
        """Create a readable color terms dictionary."""
        self.colorterms = copy.deepcopy(self.pairs)
        for second_fset in self.colorterms:
            for first_fset in self.colorterms[second_fset]:
                for filt in self.colorterms[second_fset][first_fset]:
                    print(self.colorterms[second_fset][first_fset][filt].keys())
                    self.colorterms[second_fset][first_fset][filt].pop('filter')
                    self.colorterms[second_fset][first_fset][filt].pop('colors')
                    localdic = self.colorterms[second_fset][first_fset][filt]['results']
                    for color in localdic:
                        for order in localdic[color]:
                            localdic[color][order] = np.array(localdic[color][order]['params']).tolist()

    def save_colorterms(self, output="colorterms.yaml", update=True):
        """Save results of the color term fits."""
        print("INFO: Saving results in", output)
        if os.path.exists(output) and update:
            print("INFO: Updating")
            colorterms = yaml.load(open(output, 'r'))
            colorterms.update(self.colorterms)
            results = pickle.load(open('colorfit_results.pkl', 'rb'))
            results.update(self.pairs)
        else:
            colorterms = self.colorterms
            results = self.pairs
        yaml.dump(colorterms, open(output, 'w'))
        pickle.dump(results, open('colorfit_results.pkl', 'wb'))

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
        self.init_magdiff = magdiff
        self.init_color = color
        self.magdiff = magdiff
        self.color = color
        self.kwargs = kwargs
        self.params = {}
        self.polyfits_outputs = {}

    def polyfits(self, orders="1,2,3", sigma_clip=None):
        """
        Simple polynomial fits of order 1, 2, 3.

        orders: A integer, a list of string or integer, or a string.
                e.g.: orders = 1 or "1" or "1,2,3" or ["1", "2"] or [1, 2]
        sigma_clip: Clip points above this limit and redo the fit iteratively.

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
            if order in self.polyfits_outputs:
                output = self.polyfits_outputs[order]
            else:
                output = self.polyfits_outputs[order] = {}
            # polyfit returns y = p[deg]*x**deg + p[deg-1]*x**(deg-1) + ... + p[0]
            output['params'] = list(polyfit(self.color, self.magdiff, order))
            output['ymodel'] = list(polyval(output["params"], self.color))
            output['yresiduals'] = list(self.magdiff - output['ymodel'])
            output['yresiduals_mean'] = np.mean(output['yresiduals'])
            output['yresiduals_std'] = np.std(output['yresiduals'])
            output['sigma_clip'] = np.inf if sigma_clip is None else sigma_clip
            if 'outliers' not in output:
                output['outliers'] = {'x': [], 'y': []}
            if sigma_clip is not None:
                outliers = (np.absolute(output['yresiduals']) \
                           >= sigma_clip * output['yresiduals_std'])
                while np.any(outliers):
                    # keep track of the outliers
                    output['outliers']['x'].extend(self.color[outliers])
                    output['outliers']['y'].extend(self.magdiff[outliers])
                    # and only keep good data for the next fit...
                    self.magdiff = self.magdiff[~outliers]
                    self.color = self.color[~outliers]
                    # redo the fit
                    self.polyfits(orders=str(order))
                    # check for new outliers
                    outliers = np.absolute(output['yresiduals']) \
                               >= sigma_clip * output['yresiduals_std']
                # re-compute the output using the initial input data
                output['ymodel'] = polyval(output["params"], self.init_color)
                output['yresiduals'] = self.init_magdiff - output['ymodel']
                output['yresiduals_mean'] = np.mean(output['yresiduals'])
                output['yresiduals_std'] = np.std(output['yresiduals'])
                output['sigma_clip'] = np.inf if sigma_clip is None else sigma_clip

    def plots(self, bycat_data=None, dirname="."):
        """Plot the polynomial fit results."""
        if len(self.polyfits_outputs) == 0:
            raise "INFO: You must run the polyfits method first"

        # Main figure with fits on the first axis
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111 if bycat_data is None else 121)
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
            ax.plot(self.polyfits_outputs[order]['outliers']['x'],
                    self.polyfits_outputs[order]['outliers']['y'], 'or')
        ax.legend(loc='best')

        # If by catalog data are given, add an axis where data are plotted by catalog
        if bycat_data is not None:
            ax2 = fig.add_subplot(122)
            ax2.set_xlabel(self.kwargs.get("xlabel", ""))
            ax2.set_ylabel(self.kwargs.get("ylabel", ""))
            ax2.set_title(self.kwargs.get("title", ""))
            for cdata in bycat_data:
                ax2.plot(cdata[2], cdata[0] - cdata[1], 'o', label=cdata[3])
            ax2.set_xlim(ax.get_xlim())
            ax2.set_ylim(ax.get_ylim())
            ax2.legend(loc='best')

        # Save the whole figure
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        fig.savefig("%s/%s_VS_%s.png" % \
                    (dirname,
                     self.kwargs.get("ylabel", "p1").replace(' ', ''),
                     self.kwargs.get("xlabel", "p2").replace(' ', '')))

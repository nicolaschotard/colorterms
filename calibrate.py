#!/usr/bin/env python

import os
import numpy as np

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

class TimeSpec:

    def __init__(self, day, lbda, flux, var=None):
        """ flux should be read as flux[day][lbda] """

        self.lbda = np.array(lbda)
        self.x = self.lbda
        self.flux = np.array(flux)
        self.y = self.flux
        self.var = var
        self.day = np.array(day)

    # Filters ============
    def _load_filters(self):
        """
        Load the filter names and transmission stored in the model directory.
        Also load the reference spectrum used in the model (for magnitude).
        See the following attributes:
        - FilterWheels
        - Filters
        - RefSpec
        (- VegaMags)
        """
        
        # Check first if the filters are already loaded
        if hasattr(self,'filters'):
            return

        # will load all the avalaible filters
        self.filterWheels, self.filters = load_filter_sets(self.dir, self.cards)
        
        # Load the reference spectrum:
        # Vega for SALT2 and BD17 for SALT2.2
        lbda, flux = np.loadtxt(JOIN(self.dir, self.cards['REF_SPECTRUM']),
                                unpack=True)
        self.RefSpec = OneSpec(lbda,flux)

        self.VegaMags=None
        if self.useVegaMags and self.cards.has_key('VEGAMAGS'):
            lbda,flux = np.loadtxt(JOIN(self.dir, self.cards['VEGAMAGS']),
                                   unpack=True)
            self.VegaMags=OneSpec(lbda,flux)

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
            'wrong attributes"%filt)

    # Magnitudes ============


    def mag(self, lbda, flux, var=None, step=None,
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
        elif self.VegaMags is not None:
            outmag += np.interp(filt.mean_wlength(),
                                self.VegaMags.lbda,
                                self.VegaMags.flux)

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

# ==============================================================================
# Utilities ====================================================================
# ==============================================================================

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

# SALT2 utilities====================

# Read files
def read_salt2_cards(file_name):
    """
    Read the cards SALT2 file. Return a dictionnay containing its informations
    """

    keys, values = np.loadtxt(file_name, unpack=True, dtype='string')
    cards = dict([(k[1:], v) for k, v in zip(keys, values)])

    # Now check for the version of SALT2 (GUY07 or SNLS3)
    path = file_name.replace('fitmodel.card', '')

    # Vega or BD17 spectrum
    lines = [l for l in open(path + cards['VEGA']).readlines()
             if not l.startswith('#')]
    #vk, vv = np.loadtxt(, unpack=True, dtype='string')
    #if vk[0] == '@SPECTRUM':
    if lines[0].startswith('@SPECTRUM'):
        #cards['REF_SPECTRUM'] = vv[0]
        cards['REF_SPECTRUM'] = lines[0].split(' ')[1].replace('\n', '')
    else:
        cards['REF_SPECTRUM'] = cards['VEGA']

    return cards

def load_mag_offset(path, cards):
    """
    Load the magnitude offsets to apply for the SALT2.2 mag systems.
    """
    # check the SALT2 version
    if not cards['SALT2'] == 'salt2-2-0':
        return None

    # get the magnitudes offset
    lines = np.loadtxt(JOIN(path,cards['VEGA']), dtype='string', skiprows=4)
    d=dict([(l[0],{}) for l in lines if l[0] != 'MEGACAM'])
    for l in lines:
        if l[0] == 'MEGACAM':
            continue
        d[l[0]][l[1]] = float(l[2])

    return d

def load_filter_sets(path, cards):
    """
    Load the filters include in the cards dictionnary.
    """
    FilterWheels = {}
    for syst in cards:
        if not cards[syst].startswith('Instruments') or syst == 'MEGACAMPSF':
            continue
        inst = JOIN(path, cards[syst] + '/instrument.cards')
        inst = dict(np.loadtxt(inst, dtype='string'))
        if not inst.has_key('@FILTERS'):
            continue
        FilterWheels[syst] = read_filter_wheel(JOIN(path,
                                                    JOIN(cards[syst],
                                                         inst['@FILTERS'])),
                                               cards[syst] + '/')

    # load the filter as OneSpec objects
    Filters = dict([(syst,{}) for syst in FilterWheels])
    for syst in FilterWheels:
        for filt in FilterWheels[syst]:
            # specific files for which the first row does not have the regular
            # '#' delimiter. Have to skip the first row manualy for those ones.
            # - all the MEGACAM filters
            # - the 'i_1.3' filter of the SDSS system
            # - the 'F606W' filter of the ACSWF system
            skiprows = int(syst == 'MEGACAM' or \
                           ((syst == 'SDSS') & (filt == 'i')) or \
                           ((syst == 'ACSWF') & (filt == 'F606W'))
                           )
            lbda, flux = np.loadtxt(JOIN(path, FilterWheels[syst][filt]),
                                    unpack=True, skiprows=skiprows)
            Filters[syst][filt] = OneSpec(lbda, flux)

    return FilterWheels, Filters

def read_filter_wheel(file_name, path=''):
    """
    Read the FilterWheel file to get the STANDARD set of filter names
    path: if you want to add a path before the name of the file.
    Must ends with `/`.
    """

    data = np.loadtxt(file_name, unpack=True, dtype='string')
    if len(data) == 2:
        names, files = data
    elif len(data) == 3:
        names, filts, files = data
    else:
        raise "ERROR when reading %s"%file_name

    return dict([(n, path+f) for n, f in zip(names, files)])
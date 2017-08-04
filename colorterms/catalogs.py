#!/usr/bin/env python


from glob import glob
import numpy as np
import pylab as plt
import pyfits
from pkg_resources import resource_filename
catalogs = glob(resource_filename('colorterms', 'data/catalogs/*'))
print catalogs

CATALOGS = {}


class Catalogs(object):
    pass


class Catalog(object):
    pass


def load_gunnstryker():
def load_catalogs(path="catalogs", catalog="gunnstryker", ext=".ascii",
                  desc="lbd,flux"):
    filtersets = '/'.join(path.split('/')[:-1] + ['filtersets'])
    if catalog == "gunnstryker":
        # paths
        spectra = sorted(glob(path + "/" + catalog + "/gs_*.ascii"))
        CATALOGS[catalog] = {int(sp.split('_')[1].split('.')[0]):
                             {'path': sp, 'type': 'unknown', 'name': 'unknown',
                              'lbda': None, 'flux': None, 'spec': None}
                             for sp in spectra}

        # spectrum type is any
        spectype = np.loadtxt(path + "/" + catalog + "/gsspectype.ascii",
                              dtype='str', unpack=True)
        for i, gs in enumerate(spectype[0]):
            CATALOGS[catalog][int(gs.split('_')[1].split('.')[0])]['type'] = spectype[1][i]

        # data
        for sp in spectra:
            x, y = np.loadtxt(sp, dtype='float', unpack=True)
            num = int(sp.split('_')[1].split('.')[0])
            CATALOGS[catalog][num]['lbda'] = x
            CATALOGS[catalog][num]['flux'] = y
            CATALOGS[catalog][num]['spec'] = Spectrum(x, y, fpath=filtersets)
    elif catalog == 'calspec':
        # paths
        spectra = glob(path + "/" + catalog + "/*.fits")

        # create the catalog
        CATALOGS[catalog] = {int(sp.split('_')[1].split('.')[0]):
                             {'path': sp, 'type': 'unknown', 'name': 'unknown',
                              'lbda': None, 'flux': None, 'spec': None}
                             for sp in spectra}

        # get the data
        for sp in spectra:
            spec = pyfits.open(sp)
            CATALOGS[catalog][num]['lbda'] = x
            CATALOGS[catalog][num]['flux'] = y
            CATALOGS[catalog][num]['spec'] = Spectrum(x, y, fpath=filtersets)

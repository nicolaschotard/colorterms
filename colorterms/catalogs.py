#!/usr/bin/env python


from glob import glob
import numpy as np
import pyfits
from pkg_resources import resource_filename
from . import spectools


class Catalog(object):

    def __init__(self, name, spectra, catpath=None):
        """A catalog contains a collectino of spectra.

        It is define by name."""
        self.name = name
        self.spectra = spectra
        self.catpath = catpath
        self.num_spectra = len(spectra)

        # Compute some info on the catalog
        self.mins = [min(spec.lbda) for spec in self.spectra]
        self.maxs = [max(spec.lbda) for spec in self.spectra]
        self.same_range = len(set(self.mins)) == len(set(self.maxs)) == 1
        self.min_wlength = None if not self.same_range else self.mins[0]
        self.max_wlength = None if not self.same_range else self.maxs[0]
        self.mean_wlength = None if not self.same_range else (self.mins[0] + self.maxs[0]) / 2.

    def plot_catalog(self):
        """Plot the catalog."""
        pass

    def info(self):
        print("toto")


def get_catalog_list():
    catpaths = glob(resource_filename('colorterms', 'data/catalogs/*'))
    catalogs = [catpath.split('/')[-1] for catpath in catpaths]
    return {cat: cat_path for cat, cat_path in zip(catalogs, catpaths)}


def load_catalogs():
    """Load all available catalogs and return a lost of Catalog objects."""
    catalogs = []
    catalog_names = get_catalog_list()
    for catalog in catalog_names:
        if catalog == "gunnstryker":
            catalogs.append(load_gunnstryker_catalog(catalog_names[catalog]))
        else:
            continue
    return catalogs


def load_gunnstryker_catalog(catpath):
    speclist = sorted(glob(catpath + "/gs_*.ascii"))
    spectypes = np.loadtxt(catpath + "/gsspectype.ascii", dtype='str', unpack=True)
    spectypes = {int(gs.split('_')[1].split('.')[0]): spectypes[1][i]
                 for i, gs in enumerate(spectypes[0])}
    spectra = []
    for sp in speclist:
        x, y = np.loadtxt(sp, dtype='float', unpack=True)
        name = int(sp.split('_')[1].split('.')[0])
        object_type = spectypes[name] if name in spectypes else ""
        spectra.append(spectools.Spectrum(x, y, object_name=name, object_type=object_type))
    # Order them by name/number
    sort = np.argsort([spec.object_name for spec in spectra])
    return Catalog("gunnstryker", np.array(spectra)[sort], catpath=catpath)


def load_calspec_catalog(catpath):
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
        CATALOGS[catalog][num]['spec'] = spectools.Spectrum(x, y)
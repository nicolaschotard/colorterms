#!/usr/bin/env python


from glob import glob
import numpy as np
import pylab as plt
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
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Flux')
        ax.set_title(self.name + ', %i spectra' % self.num_spectra)
        for spec in self.spectra:
            ax.plot(spec.lbda, spec.flux)
        plt.show()

    def info(self):
        print("Catalog name:", self.name)
        print("Number of spectra:", self.num_spectra)
        print("From:", self.catpath)
        print("Has same wavelength range:", self.same_range)
        if self.same_range:
            print("Wavelength info:")
            print(" - lower boundary:", self.min_wlength)
            print(" - higher boundary:", self.max_wlength)
            print(" - average value:", self.mean_wlength)


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
            print("INFO: Loading the '%s' catalog" % catalog)
            catalogs.append(load_gunnstryker_catalog(catalog_names[catalog]))
        elif catalog == "calspec":
            print("INFO: Loading the '%s' catalog" % catalog)
            catalogs.extend(load_calspec_catalog(catalog_names[catalog]))
        elif catalog == "pickles1998":
            print("INFO: Loading the '%s' catalog" % catalog)
            catalogs.extend(load_pickles1998_catalog(catalog_names[catalog]))
        else:
            print("WARNING: Loader for catalog '%s' isn't implemented yet." % catalog)
    return {cat.name: cat for cat in catalogs}


def load_gunnstryker_catalog(catpath):
    """Load the gunnstryker catalog from a catalog path and return a Catalog object."""
    # Get the list of speptra and their type
    speclist = sorted(glob(catpath + "/gs_*.ascii"))
    spectypes = np.loadtxt(catpath + "/gsspectype.ascii", dtype='str', unpack=True)
    spectypes = {int(gs.split('_')[1].split('.')[0]): spectypes[1][i]
                 for i, gs in enumerate(spectypes[0])}

    # Load each spectrum with its info
    spectra = []
    for sp in speclist:
        x, y = np.loadtxt(sp, dtype='float', unpack=True)
        object_name = int(sp.split('_')[1].split('.')[0])
        object_type = spectypes[object_name] if object_name in spectypes else ""
        spectra.append(spectools.Spectrum(x, y, object_name=object_name, object_type=object_type))

    # Order them by name/number
    sort = np.argsort([spec.object_name for spec in spectra])
    return Catalog("gunnstryker", np.array(spectra)[sort], catpath=catpath)


def load_calspec_catalog(catpath):
    """Load the calspec catalog from a catalog path and return a Catalog object.

    Five different catalogs will be return, respectively containing:
    - 'model': model spectra
    - 'stisnic': stis/nicmos data -> best quality data
    - 'fosoke': FOS+Oke data     -> second best quality data
    - 'iueoke': IUE+Oke data     - > third best quality data
    - 'bestofall': Best sptrum for each object, following the previous order
    """
    # Load the data info txt file and clean it a little
    dat = np.loadtxt(catpath + "/data_table.txt", dtype='str', delimiter='&', unpack=True)
    dat = np.array([np.array([cell.strip() for cell in col]) for col in dat])

    # Prepare for the ctalog loading
    datadict = {}
    for i, name in enumerate(dat[0]):
        # do not load the SUN spectrum
        if name == "SUN":
            continue
        bestofall = [v for v in [dat[6][i], dat[7][i], dat[8][i]] if v != ''][0]
        datadict[name] = {'name': name,
                          'type': dat[1][i],
                          'prefix': dat[4][i],
                          'model': dat[5][i],
                          'stisnic': dat[6][i],
                          'fosoke': dat[7][i],
                          'iueoke': dat[8][i],
                          'bestofall': bestofall}

    # Load catalogs
    catalogs = []
    for suffix in ["model", "stisnic", "fosoke", "iueoke", "bestofall"]:
        spectra = []
        for target in datadict:
            locd = datadict[target]  # shortcut
            if locd[suffix] == '':
                continue
            spec = pyfits.open(catpath + "/%s%s.fits" % (locd["prefix"], locd[suffix]))[1]
            spectra.append(spectools.Spectrum(spec.data.WAVELENGTH, spec.data.FLUX,
                                              object_name=locd['name'], object_type=locd['type']))
        catalogs.append(Catalog("calspec_%s" % suffix, np.array(spectra), catpath=catpath))
    return catalogs


def load_pickles1998_catalog(catpath):
    """Load the pickles1998 catalog from a catalog path and return a Catalog object."""
    # Get the list of spectrum and the corresponding type
    speclist, typelist = np.loadtxt(catpath + '/spectra_list.txt', dtype='str', unpack=True)
    # Loop over the spectra list to build the catalog
    spectra_UVILIB, spectra_UVKLIB = [], []
    for i, spec in enumerate(speclist):
        data = np.loadtxt(catpath + "/" + spec, unpack=True)  # 0: wlgth, 1: flux
        spectra = spectra_UVKLIB if spec.startswith('uk') else spectra_UVILIB
        spectra.append(spectools.Spectrum(data[0], data[1],
                                          object_name=spec.strip('.dat'),
                                          object_type=typelist[i]))
    return [Catalog("pickles1998_UVILIB", np.array(spectra_UVILIB), catpath=catpath),
            Catalog("pickles1998_UVKLIB", np.array(spectra_UVKLIB), catpath=catpath)]

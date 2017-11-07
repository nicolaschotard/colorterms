.. image:: https://landscape.io/github/nicolaschotard/colorterms/master/landscape.svg?style=flat
   :target: https://landscape.io/github/nicolaschotard/colorterms/master
   :alt: Code Health

colorterms
==========

Compute color terms that allows to go from one filter system to an other::

  colorTerms.py [options]

Available options are::

  -h, --help           show this help message and exit
  --sets SETS          Filter sets s1, s2, s3, etc. Coma separated. E.g.,
                       "sdss,megacam".Any number of filter sets can be given.
                       Fits for all combinations will be done. Set this option
                       to "all" to fit for all possible combinations between
                       available filter sets. (default: None)
  --cuts CUTS          A yaml file containing cuts to be applied on magnitudes
                       or colors.You can use the default cuts file by setting
                       this option to 'default'. (default: None)
  --sigma SIGMA        Iterative sigma clipping while fitting. (default: None)
  --catalogs CATALOGS  List of catalogs to use OR to exclude. Coma separated.
                       Add a dash at the end of your list to exclude (e.g.:
                       c1,c2,-). You can also select (or exclude) several
                       catalogs with names starting identicaly (e.g.,
                       calspec,- will exclude all calspec catalogs. (default:
                       None)
  --saveto SAVETO      Name of the file in which to save results. (default:
                       colorterms.yaml)
  --show SHOW          Show the list of available filter sets, catalogs, or
                       the content of the default yaml file containg the
                       'cuts', and exitAvailable values are: filters,
                       catalogs, cuts, and all. (default: None)

Main results are available in the `results <results>` directory.

References
==========

Here are some references about filters and colorterms previously calculated:

- `Most recent study from CADC about Megacam Vs SDSS & PanSTARRS
  filters
  <http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html>`_
  this page also contains URL to the Megacam filter bandpasses.
- `Colorterms Megacam - SDSS
  <http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/CFHTLS-SG/docs/extra/filters.html>`_. Calibration
  from `CFHT
  <http://cfht.hawaii.edu/Instruments/Imaging/MegaPrime/specsinformation.html#P2>`_
  and from `SNLS
  <http://www.astro.uvic.ca/~pritchet/SN/Calib/ColourTerms-2006Jun19/index.html#SDSScolcut>`_
  for the old `u` and `griz` bands respectively.
- `Colorterms Megacam - SDSS for the new i2 Megacam filter
  <http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/ifilt.html>`_
- `Colorterms HSC - PanSTARRS
  <https://community.lsst.org/t/pan-starrs-reference-catalog-in-lsst-format/1572>`_

Reference catalogs:

- `Gunn-Stryker <http://www.stsci.edu/hst/observatory/crds/astronomical_catalogs.html#gunn-stryker>`_

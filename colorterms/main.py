"""Main entry point for colorterms package scripts."""


import sys
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from pkg_resources import resource_filename
import yaml
from . import filtersets as Filtersets
from . import catalogs as Catalogs
from . import colorterms as Colorterms


def colorfit(argv=None):
    """Compute color terms to go from one filter set to an other."""

    description = """Compute color terms to go from one filter set to an other."""
    prog = "colorfit.py"
    usage = """%s [options]""" % prog

    parser = ArgumentParser(prog=prog, usage=usage, description=description,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--cuts', default=None,
                        help='A yaml file containing cuts to be applied on magnitudes or colors.')
    parser.add_argument('--sets', default='sdss,megacam',
                        help='Filter sets, coma separated. E.g., "')
    parser.add_argument('--show', action='store_true', default=False,
                        help="Show the list of available filter sets and catalogs along with an "
                        "exemple of what should the 'cuts' yaml file should look like, and exit")
    args = parser.parse_args(argv)

    if args.show:
        # Show the list of available filter sets
        filters = Filtersets.Filters(load=False)
        print("Available filter sets are:")
        for filt in sorted(filters.filtersets.keys()):
            print(" - %s" % filt)

        # Show the list of available catalogs
        catalogs = Catalogs.get_catalog_list()
        print("\nAvailable catalogs are:")
        for cat in sorted(catalogs):
            print(" - %s" % cat)

        # Show how the 'cuts' yaml file should look like
        cuts = open(resource_filename('colorterms', 'cuts.yaml'), 'r')
        print("\nThe 'cuts' yaml file should look like that")
        print(cuts.read())
        cuts.close()

        # Exit
        sys.exit()

    # Make sure that the input filter sets are valid
    filtersets = Filtersets.Filters(load=False)
    if len(args.sets.split(',')) != 2:
        raise IOError("You must give 2 valid sets of filters.")
    elif any([filt not in filtersets.filtersets for filt in args.sets.split(',')]):
        raise IOError("Some given filter sets do not exist (--show to get the full list).")
    else:
        filtersets = args.sets.split(',')

    # Check if a 'cuts' dictionnary has been given
    if args.cuts is not None:
        args.cuts = yaml.load(open(args.cuts, 'r'))

    # Load filter sets and catalogs
    filters = Filtersets.Filters()
    catalogs = Catalogs.load_catalogs()

    # Initialize and run the color terms computation
    colorterm = Colorterms.Colorterms(catalogs, filters)
    colorterm.compute_colorterms(filtersets[0], filtersets[1], cuts=args.cuts, sigma_clip=5)

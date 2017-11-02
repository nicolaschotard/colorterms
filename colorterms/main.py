"""Main entry point for colorterms package scripts."""


import sys
from itertools import combinations
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from pkg_resources import resource_filename
import yaml
from . import filtersets as Filtersets
from . import catalogs as Catalogs
from . import colorfits as Colorfits


def colorterms(argv=None):
    """Compute color terms to go from one filter set to an other.

    s1(f) = s2(f') + a * [s2(f') - s2(f'')] (+ other polynomial terms if asked)
    """

    description = """Compute color terms to go from one filter set to an other."""
    prog = "colorfit.py"
    usage = """%s [options]""" % prog

    parser = ArgumentParser(prog=prog, usage=usage, description=description,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sets', default=None,
                        help='Filter sets s1, s2, s3, etc. Coma separated. E.g., "sdss,megacam".'
                        'Any number of filter sets can be given. Fits for all combinations will '
                        'be done. Set this option to "all" to fit for all possible combinations '
                        'between available filter sets.')
    parser.add_argument('--cuts', default=None,
                        help='A yaml file containing cuts to be applied on magnitudes or colors.'
                        "You can use the default cuts file by setting this option to 'default'.")
    parser.add_argument('--sigma', default=None,
                        help='Iterative sigma clipping while fitting.')
    parser.add_argument('--catalogs', default=None,
                        help="List of catalogs to use OR to exclude. Coma separated."
                        " Add a dash at the end of your list to exclude (e.g.: c1,c2,-)."
                        " You can also select (or exclude) several catalogs with names "
                        "startig identicaly (e.g., calspec,- will exclude all calspec catalogs.")
    parser.add_argument('--saveto', default="colorterms.yaml",
                        help='Name of the file in which to save results.')
    parser.add_argument('--show', default=None,
                        help="Show the list of available filter sets, catalogs, or the content "
                        "of the default yaml file containg the 'cuts', and exit"
                        "Available values are: filters, catalogs, cuts, and all.")
    args = parser.parse_args(['--help']) if len(sys.argv) == 1 else parser.parse_args(argv)

    if args.show is not None:
        # Show what is asked by the user and exit
        if args.show in ("filters", "all"):
            # Show the list of available filter sets
            filters = Filtersets.Filters(load=False)
            print("Available filter sets are:")
            for filt in sorted(filters.filtersets.keys()):
                print(" - %s" % filt, filters.filtersets[filt])
            
        if args.show in ("catalogs", "all"):
            # Show the list of available catalogs
            catalogs = Catalogs.load_catalogs(verbose=False)
            print("Available catalogs are:")
            for cat in sorted(catalogs):
                print(" - %s" % cat)

        if args.show in ("cuts", "all"):
            # Show how the 'cuts' yaml file should look like
            cuts = open(resource_filename('colorterms', 'default_cuts.yaml'), 'r')
            print("Here is the default 'cuts' yaml file")
            print(cuts.read())
            cuts.close()

        if args.show not in ("filters", "catalogs", "cuts", "all"):
            print("ERROR: Available data to be shown are: 'filters', 'catalogs', 'cuts', or 'all'")

        # Exit
        sys.exit()

    # Make sure that the input filter sets are valid
    filtersets = list(Filtersets.Filters(load=False).filtersets.keys())
    if args.sets == 'all':
        # Take all filter sets
        filterset_combinations = combinations(filtersets, 2)
    elif len(args.sets.split(',')) < 2:
        raise IOError("You must give at least 2 valid sets of filters.")
    elif any([filt not in filtersets for filt in args.sets.split(',')]):
        raise IOError("Some given filter sets do not exist (get the list with `--show filters`).")
    else:
        filterset_combinations = combinations(args.sets.split(','), 2)

    # Check if a 'cuts' dictionnary has been given
    if args.cuts is not None:
        if args.cuts == "default":
            args.cuts = yaml.load(open(resource_filename('colorterms', 'default_cuts.yaml'), 'r'))
        else:
            args.cuts = yaml.load(open(args.cuts, 'r'))

    # Load filter sets and catalogs
    print("INFO: Loading filter sets")
    filters = Filtersets.Filters(verbose=False)
    print("INFO: Loading filter catalogs")
    catalogs = Catalogs.load_catalogs(verbose=False)

    # Check catalog list if given by the user
    if args.catalogs is not None:
        args.catalogs = args.catalogs.split(',')
        mode = 'exclude' if args.catalogs[-1] == '-' else 'select'

        # Expend the list of catalogs in case a short name has been given (e.g., calspec)
        args.catalogs = [cat for cat in catalogs
                        if any([cat.startswith(cc) for cc in args.catalogs])]

        # Build the final list of catalogs
        args.catalogs = [cat for cat in catalogs if cat not in args.catalogs] \
                        if mode == 'exclude' \
                           else [cat for cat in catalogs if cat not in args.catalogs]
        print("INFO: The following list of catalogs will be used:")
        for cat in sorted(args.catalogs):
            print(" -", cat)

    # Initialize and run the color terms computation for all filter set combinations
    colorterm = Colorfits.Colorterms(catalogs, filters)
    for sets in filterset_combinations:        
        colorterm.compute_colorterms(sets[0], sets[1], cuts=args.cuts, sigma_clip=float(args.sigma),
                                     verbose=False, catalogs=args.catalogs)
    colorterm.save_colorterms(args.saveto)

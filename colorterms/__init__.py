#!/usr/bin/env python

"""Initialize the colorterms package."""

import os
import glob

# Automatically import all modules (python files)
__all__ = [os.path.basename(m).replace('.py', '') for m in glob.glob("colorterms/*.py")
           if '__init__' not in m]

# Set to True if you want to import all previous modules directly
importAll = True

if importAll:
    for pkg in __all__:
        __import__(__name__ + '.' + pkg)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PyInstaller entry point for the go-enrich-tool GUI."""

import os
import sys
import multiprocessing

# Make the `src` package importable when this script is run directly from source.
_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src')
if os.path.isdir(_SRC) and _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from go_enrich.gui import main

if __name__ == "__main__":
    # Required when the bundled app uses multiprocessing (cluster --parallel).
    multiprocessing.freeze_support()
    main()

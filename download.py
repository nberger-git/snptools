#! /usr/bin/env ipython

import sys
do_snps     = True if 'snps'     in sys.argv or len(sys.argv) == 1 else False
do_genosets = True if 'genosets' in sys.argv or len(sys.argv) == 1 else False

from snptools import Downloader
Downloader().download(do_snps, do_genosets)

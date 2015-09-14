#! /usr/bin/env ipython

import sys
import snptools
snptools.merge_db(sys.argv[1:-1], sys.argv[-1])


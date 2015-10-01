#! /usr/bin/env ipython

import sys,shelve
from snptools import Config

if sys.argv[1] == 'snp' :
  filename = Config().snp_db_file
elif sys.argv[1] == 'genoset' :
  filename = Config().genoset_db_file
else :
  print "Usage : delete.py [snp|genoset] name"
  sys.exit(0)

shelf = shelve.open(filename)
if not sys.argv[2] in shelf:
  print "No entry '%s' in database" % sys.argv[2]
else :
  del shelf[sys.argv[2]]

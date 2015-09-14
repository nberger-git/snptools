#! /usr/bin/env ipython

import snptools
import sys
import os.path

logfile = 'log_genosets.txt'
if os.path.exists(logfile) :
  print 'INFO : logfile %s already exists, seems like this job has already run. Moving on.' % logfile
else:
  sys.stdout = open(logfile, 'w', 0)
  snptools.download_genosets('genosets', 'genosets.db')


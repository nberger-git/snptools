#! /usr/bin/env ipython

import snptools
import sys
import os.path

n = 500
max = int(74000/n) + 1

stdout = sys.stdout

for i in range(max - 1, 0, -1) :
  start = n*i
  stop  = n*(i+1)
  logfile = 'log_%d.txt' % i
  sys.stdout = stdout
  if os.path.exists(logfile) :
    print 'INFO : logfile %s already exists, seems like this job has already run. Moving on.' % logfile
    continue
  sys.stdout = open(logfile, 'w', 0)
  snptools.download_snp_chunk('snps', 'genotypes', start, stop, 'snp_%d.db' % i)


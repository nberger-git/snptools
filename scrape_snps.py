#! /usr/bin/env ipython

import snptools
import sys
import os.path

n = 500

if len(sys.argv) > 2 :
  min = int(sys.argv[1])
  max = int(sys.argv[2])
else :
  min = 0
  max = int(74000/n) + 1

stdout = sys.stdout

for i in range(min, max) :
  start = n*i
  stop  = n*(i+1)
  logfile = 'log_%d.txt' % i
  sys.stdout = stdout
  if os.path.exists(logfile) :
    print 'INFO : logfile %s already exists, seems like this job has already run. Moving on.' % logfile
    continue
  sys.stdout = open(logfile, 'w', 0)
  snptools.download_snp_chunk('data/snps', 'data/genotypes', start, stop, 'snp_%d.db' % i)

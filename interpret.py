#! /usr/bin/env ipython

# Usage ./interpret.py genome data/snps.db data/genosets.db output.html

import sys
import snptools
results_snps = snptools.interpret_snps(sys.argv[1], sys.argv[2])
snptools.make_snp_output(results_snps, sys.argv[2], sys.argv[4], 2)

#for r in results_snps : print r, results_snps[r]
#print len(results_snps)

results_sets = snptools.interpret_genosets(sys.argv[1], sys.argv[3], sys.argv[2])
snptools.make_genoset_output(results_sets, sys.argv[4], 2, 'a')

#print results_sets
#for r in results_sets : print r, results_sets[r]
#print len(results_sets)


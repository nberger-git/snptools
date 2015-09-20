#! /usr/bin/env ipython

# Usage ./interpret.py genome.txt output.html

import sys
from snptools import Genome, HtmlOutput
genome = Genome(sys.argv[1])

#results_snps = genome.find_snps()
#HtmlOutput().snp_page(results_snps, sys.argv[2], 2)

results_sets = genome.find_genosets()
HtmlOutput().genoset_page(results_sets, sys.argv[2], 2)

#print results_sets
#for r in results_sets : print r, results_sets[r]
#print len(results_sets)


#! /usr/bin/env ipython


import snptools
results_snps = snptools.interpret_snps('genome_test.txt', 'merge.db')
snptools.make_snp_output(results_snps, 'test.html', 2)

#for r in results_snps : print r, results_snps[r]
#print len(results_snps)

results_sets = snptools.interpret_genosets('genome_test.txt', 'genosets.db', 'merge.db')
snptools.make_genoset_output(results_sets, 'test.html', 2, 'a')

#print results_sets
#for r in results_sets : print r, results_sets[r]
#print len(results_sets)


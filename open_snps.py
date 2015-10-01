#! /usr/bin/env ipython

from snptools import Config
import shelve

config = Config()
snps = shelve.open(config.snp_db_file)

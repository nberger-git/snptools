#! /usr/bin/env python

import sys
from snptools import Data

if len(sys.argv) > 1 and sys.argv[1] == "snps" :
  Data().dump_snp_db()
elif len(sys.argv) > 1 and sys.argv[1] == "genosets" :
  Data().dump_genoset_db()
else :
  print "please specify 'snps' or 'genosets'"


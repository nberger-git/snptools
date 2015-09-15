from wikitools import wiki, category, page
import mwparserfromhell
import re
import shelve
import time
import os.path
import sys

def download_snp_list(filename = 'snps') :
  site = wiki.Wiki("http://bots.snpedia.com/api.php")                  # open snpedia
  snps = category.Category(site, "Is_a_snp")
  outfile = open(filename, 'w')
  i = 0
  for article in snps.getAllMembersGen(namespaces=[0]) :
    snp = article.title.lower()
    outfile.write(snp + '\n')
    i = i + 1
    if i % 1000 == 0 : print 'snp %5d : %20s' % (i, snp)
  outfile.close()


def download_genotype_list(filename = 'genotypes') :
  site = wiki.Wiki("http://bots.snpedia.com/api.php")                  # open snpedia
  snps = category.Category(site, "Is_a_genotype")
  outfile = open(filename, 'w')
  i = 0
  for article in snps.getAllMembersGen(namespaces=[0]) :
    genotype = article.title.lower()
    outfile.write(genotype + '\n')
    i = i + 1
    if i % 1000 == 0 : print 'genotype %5d : %20s' % (i, genotype)
  outfile.close()


def download_genoset_list(filename = 'genosets') :
  site = wiki.Wiki("http://bots.snpedia.com/api.php")                  # open snpedia
  snps = category.Category(site, "Is_a_genoset")
  outfile = open(filename, 'w')
  i = 0
  for article in snps.getAllMembersGen(namespaces=[0]) :
    genoset = article.title.lower()
    outfile.write(genoset + '\n')
    i = i + 1
    if i % 1000 == 0 : print 'genoset %5d : %20s' % (i, genoset)
  outfile.close()


def parse_name(name) :
  splits = re.split('\(|;|\)', name)
  if len(splits) > 4 : splits = splits[:-1]
  root = splits[0]
  if root[0]  == 'i'  and root[1:].isdigit() : return ['i']  + splits
  if root[:2] == 'rs' and root[2:].isdigit() : return ['rs'] + splits
  if root[:2] == 'gs' and root[2:].isdigit() : return ['gs'] + splits
  return ['other'] + splits


def get_val(t, key) :
  try:
    return str(t.get(key).value).rstrip('\n')
  except:
    print 'WARNING : could not get attribute %s from template data' % key
    #print 'INFO : could not get attribute %s from template data below \n' % key
    #print t
    return ''

def load_list(filename = 'snps', verbose = False) :
  infile = open(filename, 'r')
  #return [ name.rstrip('\n') for name in infile ]
  n_i = [0,0]
  n_rs = [0,0]
  n_gs = 0
  n_other = [0,0]
  output = []
  for name in infile :
    name = name.rstrip('\n')
    parsed = parse_name(name)
    gen = 0
    if len(parsed) > 2 : gen = 1
    if parsed[0] == 'i' :
      n_i[gen] = n_i[gen] + 1
      output.append(name)
    elif parsed[0] == 'rs' :
      n_rs[gen] = n_rs[gen] + 1
      output.append(name)
    elif parsed[0] == 'gs' :
      n_gs = n_gs + 1
      output.append(name)
    else :
      n_other[gen] = n_other[gen] + 1
      if verbose : print 'INFO: skipping ' + name
  
  print 'i-type  counts : %6d %6d' % (n_i[0], n_i[1])
  print 'rs-type counts : %6d %6d' % (n_rs[0], n_rs[1])
  print 'gs-type counts : %6d'     %  n_gs
  print 'other   counts : %6d %6d' % (n_other[0], n_other[1])
  return output


def download_snp_chunk(snpfile, genotypefile, first, last, output) :
  start_time = time.time()
  snps = load_list(snpfile, True)
  genotypes = load_list(genotypefile, True)
  num = -1
  proc = 0
  geno_dict = { geno : True for geno in genotypes }
  if os.path.exists(output) :
    print 'ERROR : DB file %s already exists, will not overwrite. Exiting now' % output
    return
  print 'INFO : writing to DB file ' + output
  shelf = shelve.open(output)
  for line in snps :
    num = num + 1
    if num < first or num >= last : continue
    name = line.rstrip('\n')
    snp = download_snp(name, geno_dict)
    if snp == None : 
      print '*** ERROR downloading SNP ' + name + ', skipping'
    elif snp != [] : 
      print 'INFO : writing SNP %s, index %d (range %d to %d)' % (name, num, first, last)
      shelf[name] = snp
      proc = proc + 1
      if proc % 10 == 0 : 
        print 'INFO : syncing DB file ' + output
        shelf.sync()
  print 'INFO : processed %d entries in %.1f seconds' % (proc, time.time() - start_time)
  print 'INFO : closing DB file ' + output
  shelf.close()
    

def download_snp(name, geno_dict = {}) :
  print 'Downloading ' + name
  parsed = parse_name(name)
  if parsed[0] != 'rs' : 
    print 'INFO : SNP %s is not rs-type, skipping' % name # only handle rsXXXXX for now
    return []
  try:
    site = wiki.Wiki("http://bots.snpedia.com/api.php")
    pagehandle = page.Page(site, name)
    snp_page = pagehandle.getWikiText()
  except:
    print 'ERROR : could not download snp %s, skipping'
    return None
  wikicode = mwparserfromhell.parse(snp_page)
  tt = wikicode.filter_templates()
  if tt[0].name != 'Rsnum\n' :
    print '*** ERROR : rs-type SNP %s has invalid record without an Rsnum entry' % name
    return None
  ok_v4 = False
  for t in tt : 
    if t.name == 'on chip ' and t.params[0] == ' 23andMe v4' :
      ok_v4 = True
      break
  if not ok_v4 :
    print 'INFO : skipping SNP %s since it is not on the 23andMe v4 chip' % name
    return []
  # Next time: add allele freq information from SNP "population diversity" template (see e.g. rs2011077)
  snp = [ #get_val(tt[0], 'rsid'), 
          get_val(tt[0], 'Chromosome'), 
          '', 
          get_val(tt[0], 'Gene'), 
          get_val(tt[0], 'position'), 
          #get_val(tt[0], 'Assembly'), 
          get_val(tt[0], 'GenomeBuild'), 
          get_val(tt[0], 'dbSNPBuild') ]
  o1 = get_val(tt[0], 'Orientation')
  o2 = get_val(tt[0], 'StabilizedOrientation')
  if len(o1) > 0 : 
    snp[1] = o1[0]
  else :
    snp[1] = '-'
  if len(o2) > 0 : 
    snp[1] = snp[1] + o2[0] 
  else :
    snp[1] = snp[1] + '-'    
  genos = []
  for i in range(1,5) :
    if tt[0].has('geno%d' % i) : 
      genotype = name + get_val(tt[0], 'geno%d' % i)
      parse_geno = parse_name(genotype)
      if len(geno_dict) and not genotype.lower() in geno_dict :
        print '*** WARNING : adding placeholder genotype for ' + genotype
        genos.append((parse_geno[2] + parse_geno[3], ['','','']))
      else :
        geno_data = download_genotype(genotype)
        genos.append((parse_geno[2] + parse_geno[3], geno_data))
  if genos == [] : return [] # no genotypes => not interesting!
  snp.append(genos)
  return snp


def download_genotype(genotype) :
  print 'Downloading ' + genotype
  site = wiki.Wiki("http://bots.snpedia.com/api.php")
  try:
    pagehandle = page.Page(site, genotype)
    snp_page = pagehandle.getWikiText()
  except:
    print '*** ERROR : genoset %s cannot be downloaded' % genotype
    return None
  wikicode = mwparserfromhell.parse(snp_page)
  tt = wikicode.filter_templates()
  repute_interp = { "good" : "g", "bad" : "b", "" : "" }
  t = tt[0]
  for t in tt :
    if t.name == 'Genotype\n' : break
  if t.name != 'Genotype\n' :
    print '*** ERROR : rs-type genotype %s has invalid record without a Genotype entry' % genotype
    return None
  repute = get_val(t, 'repute').lower()
  if repute in repute_interp : 
    repute_short = repute_interp[repute]
  else :
    repute_short = ''
  gen_data = [ get_val(t, 'magnitude'), 
               get_val(t, 'summary'),
               repute_short ]
  return gen_data


def download_genosets(genosetfile, output) :
  start_time = time.time()
  genosets = load_list(genosetfile, True)
  num = -1
  proc = 0
  if os.path.exists(output) :
    print 'ERROR : DB file %s already exists, will not overwrite. Exiting now' % output
    return
  print 'INFO : writing to DB file ' + output
  shelf = shelve.open(output)
  for line in genosets :
    num = num + 1
    name = line.rstrip('\n')
    genoset = download_genoset(name)
    if genoset == None : 
      print '*** ERROR downloading genoset ' + name + ', skipping'
    elif genoset != [] : 
      print 'INFO : writing genoset %s, index %d' % (name, num)
      shelf[name] = genoset
      proc = proc + 1
      if proc % 10 == 0 : 
        print 'INFO : syncing DB file ' + output
        shelf.sync()
  print 'INFO : processed %d entries in %.1f seconds' % (proc, time.time() - start_time)
  print 'INFO : closing DB file ' + output
  shelf.close()
    

def download_genoset(name) :
  print 'Downloading ' + name
  site = wiki.Wiki("http://bots.snpedia.com/api.php")
  try:
    pagehandle = page.Page(site, name)
    genoset_page = pagehandle.getWikiText()
    wikicode = mwparserfromhell.parse(genoset_page)
    pagehandle = page.Page(site, name + '/criteria')
    criteria_page = pagehandle.getWikiText()
  except:
    print 'ERROR : could not download genoset %s, skipping' % name
    return None
  tt = wikicode.filter_templates()
  if tt[0].name.rstrip('\n').rstrip(' ') != 'Genoset' :
    print '*** ERROR : genoset %s has invalid record without a Genoset entry' % name
    return None
  repute_interp = { "good" : "g", "bad" : "b", "" : "" }
  repute = get_val(tt[0], 'Repute').lower()
  if repute in repute_interp : 
    repute_short = repute_interp[repute]
  else :
    repute_short = ''
  genoset = [ get_val(tt[0], 'Magnitude'), 
              get_val(tt[0], 'Summary'),
              repute_short ]
  # remove # comments
  crit_text = ''
  for l in criteria_page.split('\n') : 
    print l
    toks = l.split()
    if len(toks) == 0 : continue
    if toks[0] == '#' : continue
    crit_text = crit_text + l + '\n'
  # remove HTML comments
  import lxml.html as LH
  doc = LH.fromstring(crit_text)
  criteria = doc.text_content()
  criteria = criteria.replace(' ', '').replace('\n', '').replace('\t', '')
  #import re
  #match = re.match('([a-z]*)\((.*)\)', criteria)
  #if match == None or len(match.groups()) != 2 : 
    #print '*** ERROR : genoset %s has invalid criteria %s' % (name, criteria)
    #return None
  #genoset.append(match.group(1))
  #genoset.append(match.group(2).split(','))
  genoset.append(criteria)
  return genoset


def merge_db(db_list, output) :
  print 'INFO : merging to ' + output
  output_db = shelve.open(output)
  for db_file in db_list :
    print 'INFO : ==> opening ' + db_file
    db = shelve.open(db_file, 'r')
    for snp in db : output_db[snp] = db[snp]
  output_db.close()
  

def dump_db(filename) :
  db = shelve.open(filename, 'r')
  for key in db :
    print key, db[key]

    
def read_genome(filename) :
  gen_file = open(filename, 'r')
  genome = {}
  for line in gen_file :
    tokens = line.split()
    if len(tokens) == 0 or tokens[0] == '#' : continue
    if len(tokens) != 4 or not tokens[2].isdigit() :
      print 'ERROR : invalid line "%s" while reading in genome' % line
      continue
    genome[tokens[0]] = [ tokens[3], tokens[1], tokens[2] ]
  return genome


def match_single(b1, b2, orientation = 0, rec = True) :
  u1 = b1.upper()
  u2 = b2.upper()
  if u1 == u2 :
    if orientation == -1 : 
      print 'WARNING: removing ' + u1 + u2 + ' positive match against orientation!'
    else :
      return +1
  if (u1 == 'A' and u2 == 'T') or (u1 == 'C' and u2 == 'G') :
    if orientation == +1 : 
      print 'WARNING: removing ' + u1 + u2 + ' negative match against orientation!'
    else : 
      return -1
  if u1 == 'D' and u2 == '-' : return +1
  if u1 == 'I' and re.match('[ATCG]+', u2) : return +1
  if rec :
    rev = match_single(u2, u1, orientation, False) 
    if rev != 0 : return rev
  return 0

  
def match_genotypes(t1, t2, orientation = 0) :
  if len(t1) != 1 and len(t1) != 2 : 
    print 'ERROR : invalid genotype in match :  "%s"' % t1
    return False
  if len(t2) != 1 and len(t2) != 2 : 
    print 'ERROR : invalid genotype in match "%s"' % t2
    return False
  if len(t1) != len(t2) :
    print 'ERROR : genotype length mismatch : "%s" and "%s" ' % (t1, t2)
    return False  
  if match_single(t1[0], t2[0], orientation)*match_single(t1[1], t2[1], orientation) == 1 : return True
  if match_single(t1[0], t2[1], orientation)*match_single(t1[1], t2[0], orientation) == 1 : return True
  return False


def match_genotype_list(genome_gt, ref_gt_list, orientation = 0) :
  matches = []
  for ref_gt in ref_gt_list :
    if match_genotypes(genome_gt, ref_gt, orientation) : matches.append(ref_gt)
  return matches


def interpret_snps(genome_file, db_file) :
  genome = read_genome(genome_file)
  db = shelve.open(db_file, 'r')
  results = {}
  for snpid in db :
    if not snpid in genome : 
      print 'INFO : genome SNP %s is not in the database' % snpid
      continue
    ref_snp = db[snpid]
    genome_snp = genome[snpid]
    matches = interpret_snp(snpid, ref_snp, genome_snp)
    if len(matches) == 0 : continue
    result_data = [ genome_snp[0], matches[0] ]
    ref_genotype_data = { geno[0] : geno[1] for geno in ref_snp[6] }
    print snpid
    result_data.extend(ref_genotype_data[matches[0]])
    results[snpid] = result_data
  return results
    
    
def interpret_snp(snpid, ref_snp, genome_snp) :
  if ref_snp[0] != genome_snp[1] :
    print 'ERROR : mismatch in chromosome data for snp %s : reference has "%s", genome has "%s". Skipping' % (snpid, ref_snp[0], genome_snp[1]) 
    return []
  #if ref_snp[3] != genome_snp[2] :
    #print 'WARNING : mismatch in position data for snp %s : reference has "%s", genome has "%s"' % (snpid, ref_snp[3], genome_snp[2]) 
    #continue
  ref_genotypes = [ geno[0] for geno in ref_snp[6] ]
  print 'DEBUG : testing genotype list %s for genome %s' % (str(ref_genotypes), genome_snp[0])
  if ref_snp[1] == 'p' :
    orientation = +1
  elif ref_snp[1] == 'm' :
    orientation = -1
  else:
    orientation = 0
  matches = match_genotype_list(genome_snp[0], ref_genotypes, orientation)
  print 'MATCHDEBUG', snpid, ref_genotypes, genome_snp[0], matches
  if len(matches) > 1 : 
    print 'WARNING : ambiguous matches for SNP %s : genotype %s matches %d of %s' % (snpid, genome_snp[0], len(matches), str(ref_genotypes))
    matches = []
  return matches


css_style = '''<style> 
table {
  border-collapse: collapse;
  border-top: 1px solid darkblue;
  border-bottom: 1px solid darkblue;
}
td {
font-size:1.2em;
#border: 1px solid;
border-top: 1px solid darkblue;
#border-bottom: 1px solid darkblue;
padding:2px 10px 2px 10px;
}
th {
font-size:1.4em;
text-align:left;
padding:5px 0px 5px 0px;
background-color: red;  
#color:#fff;
}
</style>
'''

def make_snp_output(snp_results, snps_file, output, mag_cutoff = 2, mode = 'w') :
  snps = shelve.open(snps_file, 'r')
  sorted_results = []
  for snpid in snp_results : 
    mag = snp_results[snpid][2]
    try:
      mag = float(mag)
    except:
      mag = 0
    sorted_results.append((snpid, mag))
  sorted_results = sorted(sorted_results, key=lambda entry: entry[1], reverse=True)
  outfile = open(output, mode)
  outfile.write(css_style)
  outfile.write('<h1> SNP Results </h1> <table>\n')
  outfile.write('<tr><td> SNP ID </td><td> Chr. </td><td> Gene </td><td> Individual Genotype </td><td> Matched Genotype </td><td> Significance </td><td> Description </td>\n')
  for snpid, mag in sorted_results :
    if mag < mag_cutoff : break # list is sorted
    results = snp_results[snpid]
    snp = snps[snpid]
    outfile.write('  <tr><td> %s </td><td> %s </td><td> %s </td> <td> %s </td>  <td> %s </td> ' % (snpid, snp[0], snp[2], results[0], results[1]))  # snpid and genotypes
    mag_style = ''
    #if results[4] == 'g' : mag_style = ' style="background-color:green'
    #elif results[4] == 'b' : mag_style = ' style="background-color:red'
    if results[4] == 'g' : mag_style = '  bgcolor="#00FF00"'
    elif results[4] == 'b' : mag_style = '  bgcolor="#FF0000"'
    outfile.write('<td%s> %s </td> <td> %s </td>' % (mag_style, results[2], results[3]))  # magnitude and summary
    outfile.write('<tr>\n')
  outfile.write('</table>\n')


def interpret_genosets(genome_file, genoset_file, snp_file) :
  genome = read_genome(genome_file)
  genosets = shelve.open(genoset_file, 'r')
  snps = shelve.open(snp_file, 'r')
  results = {}
  for genoset in genosets :
    print 'Genoset ' + genoset
    if interpret_genoset_criteria(genosets[genoset][3], genosets, snps, genome) :
      print '==> passed', genoset
      results[genoset] = genosets[genoset]
  return results

def make_genoset_output(genoset_results, output, mag_cutoff = 2, mode = 'w') :
  sorted_results = []
  for snpid in genoset_results : 
    mag = genoset_results[snpid][0]
    try:
      mag = float(mag)
    except:
      mag = 0
    sorted_results.append((snpid, mag))
  sorted_results = sorted(sorted_results, key=lambda entry: entry[1], reverse=True)
  outfile = open(output, mode)
  outfile.write(css_style)
  outfile.write('<h1> Genoset Results </h1> <table>\n')
  outfile.write('  <tr><td> Genoset ID </td><td> Significance </td><td> Description </td>\n')
  for genoid, mag in sorted_results :
    if mag < mag_cutoff : break # list is sorted
    results = genoset_results[genoid]
    outfile.write('  <tr><td> %s </td>' % genoid)
    mag_style = ''
    #if results[4] == 'g' : mag_style = ' style="background-color:green'
    #elif results[4] == 'b' : mag_style = ' style="background-color:red'
    if results[2] == 'g' : mag_style = '  bgcolor="#00FF00"'
    elif results[2] == 'b' : mag_style = '  bgcolor="#FF0000"'
    outfile.write('<td%s> %s </td> <td> %s </td>' % (mag_style, results[0], results[1]))  # magnitude and summary
    outfile.write('<tr>\n')
  outfile.write('</table>\n')
      
      
def interpret_genoset_criteria(criteria, genosets, snps, genome) :
  parsed = parse_criteria(criteria)
  print 'Parsed ' + criteria + ' ' + str(parsed)
  if parsed == None : return False
  result = evaluate_parsed_criteria(parsed, genosets, snps, genome)
  return result
  

def parse_criteria(criteria) :
  criteria = criteria.replace('\n', '').replace(' ', '').replace('\t', '')
  # criteria sometimes contains XML-style comments... remove them:
  import lxml.html as LH
  doc = LH.fromstring(criteria)
  criteria = doc.text_content()
  from pyparsing import Combine, Word, alphas, nums, Forward, Suppress, Group, ZeroOrMore, Optional
  function_call = Forward()
  identifier =  Word(alphas)
  gt = Word( alphas + '-')
  gs =  Combine('gs' + Word(nums))
  rs1 =  Combine('rs' + Word(nums) + '(' + gt + ')')
  rs2 =  Combine('rs' + Word(nums) + '(' + gt + ';' + gt + ')')
  i1 =  Combine('i' + Word(nums) + '(' + gt + ')')
  i2 =  Combine('i' + Word(nums) + '(' + gt + ';' + gt + ')')
  numb = Word(nums)
  arg = function_call | rs1 | rs2 | i1 | i2 | gs | numb
  function_call <<= Group(identifier + Suppress('(') + Group(ZeroOrMore(arg + Suppress(",")) + Optional(arg)) + Suppress(')'))
  try:
    result = arg.parseString(criteria).asList()[0]
  except :
    print 'ERROR : could not parse criteria %s' % criteria
    return None
  return result
  #return function_call.parseString(criteria).asList()[0]


def evaluate_parsed_criteria(parsed, genosets, snps, genome) :
  print 'EVP : ', parsed
  if type(parsed) == str : # simple case
    import re
    rs1 = re.match('rs([0-9]*)\(([A-Z,-]*)\)', parsed)
    if rs1 != None :
      rsid = 'rs' + rs1.group(1)
      if not rsid in genome : return False
      genome_snp = genome[rsid]
      if len(genome_snp[0]) == 1 : genome_snp[0] = genome_snp[0] + ' ' # can occur for male X-chromosome genotypes
      if not rsid in snps : return False
      snpdb_snp = snps[rsid]
      if snpdb_snp[1] == 'p' :
        orientation = +1
      elif snpdb_snp[1] == 'm' :
        orientation = -1
      else:
        orientation = 0
      req = rs1.group(2)
      print 'compare', rsid, genome_snp[0], req, orientation
      print '->', match_single(genome_snp[0][0], req, orientation) != 0 or match_single(genome_snp[0][1], req, orientation) != 0
      return match_single(genome_snp[0][0], req, orientation) != 0 or match_single(genome_snp[0][1], req, orientation) != 0
    rs2 = re.match('rs([0-9]*)\(([A-Z,-]*);([A-Z,-]*)\)', parsed)
    if rs2 != None :
      rsid = 'rs' + rs2.group(1)
      if not rsid in genome : return False
      genome_snp = genome[rsid]
      if len(genome_snp[0]) == 1 : genome_snp[0] = genome_snp[0] + ' ' # can occur for male X-chromosome genotypes
      if not rsid in snps : return False
      snpdb_snp = snps[rsid]
      if snpdb_snp[1] == 'p' :
        orientation = +1
      elif snpdb_snp[1] == 'm' :
        orientation = -1
      else:
        orientation = 0
      req = rs2.group(2) + rs2.group(3)
      print 'compare', genome_snp[0], req, orientation
      print '->', match_genotypes(genome_snp[0], req, orientation)
      return match_genotypes(genome_snp[0], req, orientation)
    gs  = re.match('gs([0-9]*)', parsed)
    if gs != None : 
      if not parsed in genosets : return False
      return interpret_genoset_criteria(genosets[parsed][3], genosets, snps, genome)
    print 'ERROR : unsupported string pattern %s in genoset evaluation' % parsed
  # now recursive case
  if len(parsed) != 2 :
    print 'ERROR : parsed criteria %s do not have 2 items (function and argument list) as expected' % str(parsed)
    return None
  if parsed[0] == 'atleast' :
    if not parsed[1][0].isdigit :
      print 'ERROR : atleast arguments do not start with a number as expected, got %2' % parsed[1][0]
      return None
    minval = int(parsed[1][0])
    results = [ evaluate_parsed_criteria(item, genosets, snps, genome) for item in parsed[1][1:] ]
    print results
    print '->', sum(1 for item in results if item), '>=?', minval, '==>', (sum(1 for item in results if item) >= minval)
    return (sum(1 for item in results if item) >= minval)
  results = [ evaluate_parsed_criteria(item, genosets, snps, genome) for item in parsed[1] ] 
  print 'EVP : results = ', results, "'" + parsed[0] + "'"
  if parsed[0] == 'and' : print '->', all(results)
  if parsed[0] == 'or'  : print '->', any(results)
  if parsed[0] == 'and' : return all(results)
  if parsed[0] == 'or'  : return any(results)
  if parsed[0] == 'not' :
    if len(results) != 1 :
      print 'ERROR : passing %d arguments to "not" in %s, should be one' % (len(parsed[1]), str(parsed))
      return None
    print '->', not results[0]
    return not results[0]
  print 'ERROR : unknown operator %s' % parsed[0]
  return None


def make_snp_db_output(snp_db_file, output, verbose = False, mode = 'w') :
  db = shelve.open(snp_db_file, 'r')
  outfile = open(output, mode)
  outfile.write(css_style)
  #rs11155133 ['6', 'p', 'LOC102723724', '140848688', '38.1', '141', [('AA', ['0', 'common in complete genomics', 'g'])]]
  outfile.write('<h1> SNP Database </h1> <table>\n')
  if verbose :
    outfile.write('<tr><td> SNP </td> <td> Chromosome </td> <td> Orientation </td>  <td> Gene </td> <td> Location </td> <td> Genome Build </td>  <td> dbSNP Build </td> <td> Genotype </td> <td> Significance </td>  <td> Description </td></tr>\n')
  else:
    outfile.write('<tr><td> SNP </td> <td> Chromosome </td>  <td> Gene </td> <td> Genotype </td> <td> Significance </td>  <td> Description </td>\n')
  for snpid in db :
    results = db[snpid]
    td = '<td rowspan="%d">' % len(results[6])
    if verbose :
      outfile.write('<tr>%s %s </td> %s %s </td>  %s %s </td> %s %s </td> %s %s </td>  %s %s </td> %s %s </td>' % 
                    (td, snpid, td, results[0], td, results[1], td, results[2], td, results[3], td, results[4], td, results[5]))
    else:
      outfile.write('<tr> %s %s </td> %s %s </td>  %s %s </td>' % 
                    (td, snpid, td, results[0], td, results[2]))
    for gt in results[6] :
      mag_style = ''
      if gt[1][2] == 'g' : mag_style = ' bgcolor="#00FF00"'
      elif gt[1][2] == 'b' : mag_style = ' bgcolor="#FF0000"'
      outfile.write('<td> %s </td> <td%s> %s </td>  <td> %s </td></tr>\n' % (gt[0], mag_style, gt[1][0], gt[1][1]))
  outfile.write('</table>\n')


def make_genoset_db_output(genoset_db_file, output, mode = 'w') :
  db = shelve.open(genoset_db_file, 'r')
  outfile = open(output, mode)
  outfile.write(css_style)
  #gs269 ['2', 'APOE E2/E3', 'g', 'and(rs429358(T;T),rs7412(C;T))']
  outfile.write('<h1> Genoset Database </h1> <table>\n')
  outfile.write('<tr><td> Genoset </td> <td> Significance </td> <td> Description </td> <td> Criteria </td></tr>\n')
  for geno in db :
    results = db[geno]
    mag_style = ''
    if results[2] == 'g' : mag_style = ' bgcolor="#00FF00"'
    elif results[2] == 'b' : mag_style = ' bgcolor="#FF0000"'
    outfile.write('<td> %s </td> <td%s> %s </td>  <td> %s </td><td> %s </td></tr>\n' % (geno, mag_style, results[0], results[1], results[3]))
  outfile.write('</table>\n')

from wikitools import wiki, category, page, NoPage
import mwparserfromhell
import re
import shelve
import time
import os
import sys
import ConfigParser


config_locations = [ os.path.join(loc, 'snptools.conf') for loc in [ os.curdir, os.path.expanduser("~/.snptools") ] ]


class Config :  
  def __init__(self, locs = config_locations) :
    config = ConfigParser.SafeConfigParser()
    found = config.read(locs)
    if found == [] :
      print 'ERROR : could not locate config file, tried ' + ', '.join(locs)
      sys.exit(1)
    try :
      self.snp_db_file = os.path.expanduser(config.get('snps', 'db_file'))
      self.snp_list_file = os.path.expanduser(config.get('snps', 'list_file'))
      self.genotype_list_file = os.path.expanduser(config.get('genotypes', 'list_file'))
      self.genoset_list_file = os.path.expanduser(config.get('genosets', 'list_file'))
      self.genoset_db_file = os.path.expanduser(config.get('genosets', 'db_file'))
      self.site = config.get('wiki', 'site')
      self.use_get_requests = config.get('wiki', 'use_get_requests') in ['True', 'true', 'yes', '1' ]
      self.css_style = config.get('html', 'css_style')
      self.verbose = config.get('output', 'verbose')
    except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as error :
      error.message = error.message + '\nERROR : cannot read config file %s, exiting' % found[0]
      raise
    message = '%s file %s is not found, please run the data fetching utility first'


class Data :
  def __init__(self, config = None) :
    self.config = config if config != None else Config()
    self.l_snps = None
    self.l_genosets = None
    self.l_snp_list = None
    self.l_genotype_list = None
    self.l_genoset_list = None

  def snps(self) :
    if not os.path.exists(self.config.snp_db_file) : raise IOError, message % ('SNP database', self.config.snp_db_file)
    if self.l_snps == None : self.l_snps = shelve.open(self.config.snp_db_file, 'r')
    return self.l_snps
  
  def genosets(self) :
    if not os.path.exists(self.config.genoset_db_file) : raise IOError, message % ('Genoset database', self.config.genoset_db_file)
    if self.l_genosets == None : self.l_genosets = shelve.open(self.config.genoset_db_file, 'r')
    return self.l_genosets

  def snp(self, snp_name) :
    try :
      return SNP(snp_name).deserialize(self.snps()[snp_name])
    except KeyError :
      return None

  def genoset(self, genoset_name) :
    try :
      return SNPGenoset(name=genoset_name, data=self).deserialize(self.genosets()[genoset_name])
    except KeyError :
      return None

  def load_list(self, filename) :
    try:
      infile = open(filename, 'r')
    except IOError as error :
      print error
      return None
    output = []
    for name in infile : output.append(name.rstrip('\n'))
    return output
    
  def snp_list(self) :
    if not os.path.exists(self.config.snp_list_file) : raise IOError, message % ('SNP list', self.config.snp_list_file)
    if self.l_snp_list == None : self.l_snp_list = self.load_list(self.config.snp_list_file)
    return self.l_snp_list

  def genotype_list(self) :
    if not os.path.exists(self.config.genotype_list_file) : raise IOError, message % ('Genotype list', self.config.genotype_list_file)
    if self.l_genotype_list == None : self.l_genotype_list = self.load_list(self.config.genotype_list_file)
    return self.l_genotype_list

  def genoset_list(self) :
    if not os.path.exists(self.config.genoset_list_file) : raise IOError, message % ('Genoset list', self.config.genoset_list_file)
    if self.l_genoset_list == None : self.l_genoset_list = self.load_list(self.config.genoset_list_file)
    return self.l_genoset_list

  def dump_snp_db(self) : self.dump_db(self.snps())
  def dump_genoset_db(self) : self.dump_db(self.genosets())
  def dump_db(self, db) :
    for key in db :
      print key, db[key]


class Wiki :
  def __init__(self, config = None) :
    self.config = config if config != None else Config()
    print 'INFO: opening wiki site %s' % self.config.site
    self.site = wiki.Wiki(self.config.site)

  def get_category(self, catname) :
    cat = category.Category(self.site, catname)
    items = []
    for article in cat.getAllMembersGen(namespaces=[0]) :
      items.append(article.title.lower())
      if len(items) % 1000 == 0 : print 'Downloading item %5d : %20s' % (len(items), items[-1])
    return items
  
  def get_text(self, name) :
    try:
      capname = name[0].upper() + name[1:]
      pagehandle = page.Page(self.site, capname, False, False)
      try:
        return pagehandle.getWikiText(getrequest=self.config.use_get_requests)
      except :
        # For unpatched wikitools versions
        return pagehandle.getWikiText()
    except NoPage :
      print 'ERROR : page %s is not found!' % capname
      return None

  def get_templates(self, name) :
    text = self.get_text(name)
    if not text : return []
    wikicode = mwparserfromhell.parse(text)
    return wikicode.filter_templates()

  def get_template(self, templates, template_name) :
    for t in templates :
      if t.name == template_name : return t
    return None
  
  def get_val(self, templates, template_name, key, warn = True) :
    t = self.get_template(templates, template_name)
    if t == None :
      print 'WARNING : no template found with name %s' % template_name
      return None
    try :
      return str(t.get(key).value).rstrip('\n')
    except:
      if warn : print 'WARNING : could not get attribute %s from template %s' % (key, template_name)
      return None
  

class Downloader :
  def __init__(self, config = None) :
    self.config = config if config != None else Config()
    self.wiki = Wiki(config)
    self.data = Data()
    
  def download_category(self, catname, output_file, overwrite = False) :
    if os.path.exists(output_file) :
      print 'INFO : output file %s already exists, skipping' % output_file
      return True
    print 'INFO : downloading category %s to file %s.' % (catname, output_file)
    results = [ item for item in self.wiki.get_category(catname) ]
    if os.path.exists(output_file) :
      print 'INFO : output file %s was created while data was retrieved, skipping' % output_file
      return True
    outfile = open(output_file, 'w')
    for item in results : outfile.write(item + '\n')
    outfile.close()
    return True
  
  def download_snp(self, snp_name, geno_dict = {}) :
    snp = SNP(snp_name)
    templates = self.wiki.get_templates(snp_name)
    if self.wiki.get_template(templates, 'Rsnum\n') :
      snp.chromosome = self.wiki.get_val(templates, 'Rsnum\n', 'Chromosome')
      snp.gene = self.wiki.get_val(templates, 'Rsnum\n', 'Gene')
      position = self.wiki.get_val(templates, 'Rsnum\n', 'position')
      snp.position = int(position) if position and position.isdigit() else None
      snp.genome_build = self.wiki.get_val(templates, 'Rsnum\n', 'GenomeBuild')
      snp.dbSNP_build = self.wiki.get_val(templates, 'Rsnum\n', 'dbSNPBuild')
      snp.orientation = self.wiki.get_val(templates, 'Rsnum\n', 'Orientation')
      snp.stabilized_orientation = self.wiki.get_val(templates, 'Rsnum\n', 'StabilizedOrientation')
    elif self.wiki.get_template(templates, '23andMe SNP\n') :
      snp.chromosome = self.wiki.get_val(templates, '23andMe SNP\n', 'Chromosome')
      snp.gene = self.wiki.get_val(templates, '23andMe SNP\n', 'Gene')
      position = self.wiki.get_val(templates, '23andMe SNP\n', 'position')
      snp.position = int(position) if position and position.isdigit() else None
      snp.genome_build = self.wiki.get_val(templates, '23andMe SNP\n', 'GenomeBuild')
      snp.dbSNP_build = None
      snp.orientation = None
      snp.stabilized_orientation = None
    snp.genotypes = {}
    for i in range(1,4) :
      val = self.wiki.get_val(templates, 'Rsnum\n', 'geno%d' % i, False) 
      if val :
        genotype_name = snp_name + val
        if len(geno_dict) and not genotype_name.lower() in geno_dict :
          print '*** WARNING : genotype %s not in list, adding placeholder instead ' % genotype_name
          genotype = SNPGenotype(snp, val)
        else :
          genotype = self.download_genotype(snp, val)
        snp.genotypes[val] = genotype
    popd = self.get_populations(templates)
    for geno in popd :
      if geno in snp.genotypes :
        snp.genotypes[geno].frequency = popd[geno]
    if snp.genotypes == {} and self.wiki.get_template(templates, 'hgsnp\n') :
      anc_all = self.wiki.get_val(templates, 'hgsnp\n', 'ancestral_allele')
      der_all = self.wiki.get_val(templates, 'hgsnp\n', 'derived_allele')
      try :
        summary = ','.join([ s.rstrip('\n') for s in self.wiki.get_template(templates, 'hgsnp\n').split("|")[1:7] ])
      except:
        summary = ''
      if anc_all : 
        anc_geno = SNPGenotype(snp, anc_all)
        anc_geno.summary = summary
        snp.genotypes[anc_all] = anc_geno
      if der_all : 
        der_geno = SNPGenotype(snp, der_all)
        der_geno.summary = summary
        snp.genotypes[der_all] = der_geno
    return snp

  def download_genotype(self, snp, genotype_name) :
    genotype = SNPGenotype(snp, genotype_name)
    templates = self.wiki.get_templates(genotype.full_name())
    genotype.magnitude = self.wiki.get_val(templates, 'Genotype\n', 'magnitude')
    try :
      genotype.magnitude = float(genotype.magnitude)
    except :
      genotype.magnitude = None
    genotype.summary   = self.wiki.get_val(templates, 'Genotype\n', 'summary')
    genotype.repute = self.wiki.get_val(templates, 'Genotype\n', 'repute').lower()
    if genotype.magnitude != 'good' and genotype.magnitude != 'bad' : genotype.magnitude = None
    return genotype

  def download_genoset(self, genoset_name) :
    genoset = SNPGenoset(genoset_name)
    templates = self.wiki.get_templates(genoset_name)
    genoset.repute = self.wikiget_val(templates, 'Genoset', 'Repute').lower()
    if genoset.repute != 'good' and genoset.repute != 'bad' : genoset.repute = None
    genoset.magnitude = get_val(templates, 'Genoset', 'Magnitude')
    try : 
      genoset.magnitude = int(genoset.magnitude) 
    except : 
      genoset.magnitude = None
    genoset.summary   = get_val(templates, 'Genoset', 'Summary')
    criteria_page = self.wiki.get_text(genoset_name + '/criteria')
    # remove # comments
    criteria_text = ''
    for l in criteria_page.split('\n') : 
      toks = l.split()
      if len(toks) == 0 : continue
      if toks[0] == '#' : continue
      criteria_text = criteria_text + l + '\n'
    # remove HTML comments
    import lxml.html as LH
    doc = LH.fromstring(criteria_text)
    criteria = doc.text_content()
    criteria = criteria.replace(' ', '').replace('\n', '').replace('\t', '')
    genoset.criteria = append(criteria)
    return genoset

  def get_populations(self, templates) :
    pops = []
    for t in templates :
      if t.name != ' population diversity\n' : continue
      for line in t.split('\n') :
        tokens = [ s.strip(' ') for s in line.split('|') ]
        if tokens[1] != 'CEU' : continue
        pops = tokens[2:]
        break
    popd = {}
    try :
      for i, pop in enumerate(pops) :
        popd[t.get('geno%d' % (i+1)).value.rstrip('\n')] = float(p)
    except:
      print 'ERROR : could not retrieve population diversity values'
    return popd
    
  def download_snps(self) :
    start_time = time.time()
    snps = self.data.snp_list()
    genotypes = self.data.genotype_list()
    geno_dict = { geno : True for geno in genotypes }
    processed = 0
    print 'INFO : writing to DB file ' + self.config.snp_db_file
    shelf = shelve.open(self.config.snp_db_file)
    interrupt = False
    try :
      for num, line in enumerate(snps) :
        name = line.rstrip('\n')
        if name in shelf : continue
        print 'INFO : downloading %s' % name
        snp = self.download_snp(name, geno_dict)
        if snp == None : 
          print 'ERROR downloading SNP ' + name + ', skipping for now'
          continue
        print 'INFO : writing SNP %s, index %d of %d' % (name, num, len(snps))
        shelf[name] = snp.serialize()
        processed = processed + 1
    except KeyboardInterrupt :
      print 'INFO: SNP download was interrupted cleanly'
      interrupt = True
    print 'INFO : processed %d entries in %.1f seconds' % (processed, time.time() - start_time)
    print 'INFO : closing DB file ' + self.config.snp_db_file
    shelf.close()
    return not interrupt

  def download_genosets(self) :
    start_time = time.time()
    genosets = self.data.genoset_list()
    processed = 0
    print 'INFO : writing to DB file ' + self.config.genoset_db_file
    shelf = shelve.open(self.config.genoset_db_file)
    interrupt = False
    try :
      for num, line in enumerate(genosets) :
        name = line.rstrip('\n')
        if name in shelf : continue
        print 'INFO : downloading %s' % name
        genoset = self.download_genoset(name)
        if genoset == None : 
          print 'ERROR downloading genoset ' + name + ', skipping for now'
          continue
        print 'INFO : writing genoset %s, index %d of %d' % (name, num, len(genosets))
        shelf[name] = genoset.serialize()
        processed = processed + 1
    except KeyboardInterrupt :
      print 'INFO: genoset download was interrupted cleanly'
      interrupt = True
    print 'INFO : processed %d entries in %.1f seconds' % (processed, time.time() - start_time)
    print 'INFO : closing DB file ' + self.config.genoset_db_file
    shelf.close()
    return not interrupt

  def download(self) :
    try :
      if self.download_category('Is_a_snp',      self.config.snp_list_file) and \
         self.download_category('Is_a_genotype', self.config.genotype_list_file) and \
         self.download_category('Is_a_genoset',  self.config.genoset_list_file) and \
         self.download_snps() and \
         self.download_genosets() : return True
    except KeyboardInterrupt :
      pass
    print 'INFO : download process was interrupted cleanly, leaving now'
    
      
class Allele :
  def __init__(self, allele, orientation) :
    self.orientation = orientation if orientation == 'plus' or orientation == 'minus' else None
    upal = allele.upper()
    if re.match('[ATCG]+', upal) : self.allele = upal # for I/D genotypes, the 'I' case may have mroe than 1 base
    elif upal == 'I' : self.allele = '.'
    elif upal == '-' : self.allele = '-'
    elif upal == 'D' : self.allele = '-'
    else :
      print 'ERROR : invalid allele %s in genotype' % allele
      self.allele = None
  
  def plus_allele(self) :
    if self.orientation == None : return None
    if self.orientation == 'plus' : return self.allele
    if self.allele == '.' : return '.'
    if self.allele == '-' : return '-'
    flipmap = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C' }
    return ''.join(flipmap[c] for c in self.allele)
  
  def match(self, other) :
    u1 = self.plus_allele()
    u2 = other.plus_allele()
    print 'yyy comparing', u1, u2
    if u1 == u2 : return True
    if u1 == '.' and u2 != '-' : return True
    if u2 == '.' and u1 != '-' : return True
    return False

  @staticmethod
  def split(genotype, orientation) :
    if len(genotype) <= 2 : return [ Allele(c, orientation) for c in genotype ] # A or GG
    if genotype[0] == '-' : return [ Allele(genotype[0], orientation), Allele(genotype[1:], orientation) ] # -ACG
    l = len(genotype)
    if l % 2 == 0 : return [Allele(genotype[0:l/2], orientation), Allele(genotype[l/2:], orientation) ] # ACGACG
    return None


class SNP :
  def __init__(self, name) :
    self.name = name
    self.chromosome = None
    self.gene = None
    self.position = None
    self.genome_build = None
    self.dbSNP_build = None
    self.orientation = None
    self.stabilized_orientation = None
    self.frequency = None
    self.genotypes = {}
    
  def match(self, genotype) :
    print 'xxx comparing', genotype.genotype, 'to snp', self.name
    if genotype.chromosome != None and self.chromosome != None and genotype.chromosome != self.chromosome : return []
    if genotype.position != None and self.position != None and genotype.position != self.position : return []
    matches = [ self.genotypes[g] for g in self.genotypes if self.genotypes[g].match(genotype) ]
    print 'xxx matches:', ','.join(g.genotype for g in matches)
    return matches

  def serialize(self) :
    record = [ self.chromosome, self.gene, self.position, self.genome_build, self.dbSNP_build ]
    ori = '-'
    sor = '-'
    if self.orientation and len(self.orientation) > 0 : ori = self.orientation[0]
    if self.stabilized_orientation and len(self.stabilized_orientation) > 0 : sor = self.stabilized_orientation[0]
    record.append(ori + sor)
    genotype_records = {}
    for geno in self.genotypes :
      genotype_records[geno] = self.genotypes[geno].serialize()
    record.append(genotype_records)
    return record
  
  def deserialize(self, record) :
    # HACK for now
    self.chromosome = record[0]
    self.gene = record[2]
    self.position = int(record[3]) if record[3].isdigit() else None
    self.genome_build = record[4]
    self.dbSNP_build = record[5]
    ori = record[1][0] if len(record[1]) > 0 else None
    if ori == 'p' : self.orientation = 'plus'
    elif ori == 'm' : self.orientation = 'minus'
    self.stabilized_orientation = self.orientation
    for geno in record[6] :
      self.genotypes[geno[0]] = SNPGenotype(self, geno[0]).deserialize(geno[1])
    return self
    #end HACK
    try :
      self.chromosome = record[0]
      self.gene = record[1]
      self.position = record[2]
      self.genome_build = record[3]
      self.dbSNP_build = record[4]
      ori = record[5][0]
      if ori == 'p' : self.orientation = 'plus'
      elif ori == 'm' : self.orientation = 'minus'
      else : self.orientation = None
      sor = record[5][1]
      if sor == 'p' : self.stabilized_orientation = 'plus'
      elif sor == 'm' : self.stabilized_orientation = 'minus'
      else : self.stabilized_orientation = None
      for geno in record[6] :
        self.genotypes[geno] = SNPGenotype(self, geno).deserialize(record[6][geno])
    except KeyError as error :
      error.message = error.message + '\nERROR : could not deserialize snp record %s' % record
      raise
    return self
      
class SNPGenotype() :
  def __init__(self, snp, genotype) :
    if type(genotype) == list :
      self.alleles = genotype
    else :
      self.alleles = self.parse_alleles(genotype)
    self.snp = snp
    self.magnitude = None
    self.summary = None
    self.repute = None
        
  def parse_alleles(self, raw_genotype) :
    match = re.match('\((([A-Z]|-)*)(;(([A-Z]|-)*))?\)', raw_genotype)
    if match : return [ match.group(1) ].append(match.group(6) if match.group(6) else [])
    match = re.match('^([A-Z]|-)([A-Z]|-)$', raw_genotype)
    if match : return [ match.group(0), match.group(1) ]
    match = re.match('^([A-Z]|-)$', raw_genotype)
    if match : return [ match.group(0) ]
    return None

  def full_name(self) : 
    return self.snp.name + self.genotype

  def alleles(self) :
    return [ Allele(a, snp.orientation) for a in self.alleles ]

  def match(self, genotype) :
    ref_alleles = self.alleles()
    gen_alleles = genotype.alleles()
    if len(ref_alleles) == 1 : # snp_genotypes specifies 1 allele : can match either of ours
      for a in gen_alleles : 
        if a.match(ref_alleles[0]) : return True
      return False
    if len(ref_alleles) == 2 : # snp_genotypes specifies 2 allele : must match both of ours
      if len(gen_alleles) != 2 : return False
      print 'xxx compare genotype gen', ','.join(al.allele for al in gen_alleles)
      print 'xxx compare genotype snp', ','.join(al.allele for al in ref_alleles)
      if gen_alleles[0].match(ref_alleles[0]) and gen_alleles[1].match(ref_alleles[1]) : return True
      if gen_alleles[0].match(ref_alleles[1]) and gen_alleles[1].match(ref_alleles[0]) : return True
      return False
    return False

  def match_genome(self, genome) :
    genotype = genome.genotype(self.snp.name)
    return self.match(genotype) if genotype else False

  def serialize(self) :
    if self.repute == 'good' : rep = True 
    elif self.repute == 'bad' : rep = False
    else : rep = None
    return [ self.magnitude, rep, self.summary ]
  
  def deserialize(self, record) :
    # HACK for now
    try:
      self.magnitude = float(record[0])
    except:
      self.magnitude = None
    self.summary = record[1]
    rep = record[2]
    if rep == 'g' :  self.repute = 'good'
    elif rep == 'b' : self.repute = 'bad'
    else : self.repute = None
    return self
    #end HACK
    try :
      try:
        self.magnitude = float(record[0])
      except : self.magnitude = None
      if record[1] == True : self.repute = 'good'
      elif record[1] == False : self.repute = 'bad'
      else : self.repute = None
      self.summary = record[2]
    except KeyError as error :
      error.message = error.message + '\nERROR : could not deserialize genotype record %s' % record
      raise
    return self


class SNPGenoset() :
  def __init__(self, name = None, criteria = None, parsed = None, data = None) :
    self.name = name
    self.magnitude = None
    self.repute = None
    self.criteria = criteria
    self.summary = None
    self.parsed = parsed
    self.data = data if data != None else Data()

  def parsed_criteria(self) :
    print 'qqq1', self.name, self.criteria, self.parsed
    if self.parsed != None : return self.parsed
    from pyparsing import Combine, Word, alphas, nums, Forward, Suppress, Group, ZeroOrMore, Optional, ParseException
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
      result = arg.parseString(self.criteria).asList()[0]
    except ParseException as error :
      #error.message = error.message + '\nERROR : could not parse criteria %s' % self.criteria
      #raise
      print  'ERROR : could not parse criteria %s' % self.criteria
      result = []
    self.parsed = result
    print 'qqq', result
    return self.parsed

  def parse_snp_genotype(self) : 
    match = re.match('(rs|i)([0-9]*)\((([A-Z]|-)*)(;(([A-Z]|-)*))?\)', self.parsed_criteria())
    if match == None : return None
    snp = self.data.snp(match.group(1) + match.group(2))
    if snp == None : 
      print 'WARNING : while evaluating criteria %s, encountered unknown SNP %s' % (self.parsed_criteria(), match.group(1) + match.group(2))
      return False
    print 'ggg', snp.name, match.groups()
    return SNPGenotype(snp, match.group(3) + (match.group(6) if match.group(6) else ''))

  def parse_genoset(self) :
    gs  = re.match('gs([0-9]*)', self.parsed_criteria())
    if gs == None : return None
    return self.data.genoset(self.parsed_criteria())
    
  def match_genome(self, genome) :
    if type(self.parsed_criteria()) == str : # simple case : a single expression
      snp_genotype = self.parse_snp_genotype()
      if snp_genotype == False : return False # SNP not present, skip
      if snp_genotype != None : return snp_genotype.match_genome(genome)
      # it is None => could not parse
      genoset = self.parse_genoset()
      if genoset : return genoset.match_genome(genome)
      print 'ERROR : unsupported string pattern %s in genoset evaluation' % self.parsed_criteria()
      return False
    # If not, we have a function
    if len(self.parsed_criteria()) != 2 :
      print 'ERROR : parsed criteria %s do not have 2 items (function and argument list) as expected' % str(self.parsed_criteria())
      return False
    func = self.parsed_criteria()[0]
    args = self.parsed_criteria()[1]
    if func == 'atleast' :
      if not args[0].isdigit :
        print 'ERROR : atleast arguments do not start with a number as expected, got %2' % args[0]
        return False
      minval = int(args[0])
      results = [ SNPGenoset(parsed=item, data=self.data).match_genome(genome) for item in args[1:] ]
      print results
      print '->', sum(1 for item in results if item), '>=?', minval, '==>', (sum(1 for item in results if item) >= minval)
      return (sum(1 for item in results if item) >= minval)
    # other functions now:
    print '???', args
    results = [ SNPGenoset(parsed=item, data=self.data).match_genome(genome) for item in args ] 
    print 'EVP : results = ', results, "'" + func + "'"
    if func == 'and' : print '->', all(results)
    if func == 'or'  : print '->', any(results)
    if func == 'and' : return all(results)
    if func == 'or'  : return any(results)
    if func == 'not' :
      if len(results) != 1 :
        print 'ERROR : passing %d arguments to "not" in %s, should be one' % (len(args), str(self.parsed_criteria()))
        return None
      print '->', not results[0]
      return not results[0]
    print 'ERROR : unknown operator %s' % func
    return None

  def serialize(self) :
    if self.repute == 'good' : rep = True 
    elif repute == 'bad' : rep = False
    else : rep = None
    return [ self.magnitude, rep, self.criteria, self.summary ]
  
  def deserialize(self, record) :
    # HACK for now
    print 'rrr', self.name, record
    try:
      self.magnitude = float(record[0])
    except:
      self.magnitude = None
    if record[2] == 'g' : self.repute = 'good'
    elif record[2] == 'b' : self.repute = 'bad'
    else : self.repute = None
    self.criteria = record[3]
    self.summary = record[1]    
    return self
    #end HACK
    try :
      self.magnitude = record[0]
      if record[1] == True : self.repute = 'good'
      elif record[1] == False : self.repute = 'bad'
      else : self.repute = None
      self.criteria = record[2]
      self.summary = record[3]
    except KeyError as error :
      error.message = error.message + '\nERROR : could not deserialize genoset record %s' % record
      raise
    return self
  

class Genome :
  def __init__(self, filename, data = None) :
    self.data = data if data != None else Data()
    self.genotypes = {}
    gen_file = open(filename, 'r')
    for line in gen_file :
      tokens = line.split()
      if len(tokens) == 0 or tokens[0] == '#' : continue
      if len(tokens) != 4 :
        print '*** ERROR : invalid line "%s" while reading in genome' % line
        continue
      # 23andme format:
      # rsid  chromosome      position        genotype
      # rs4477212       1       82154   AA
      # orientation = plus since 23andm3 always reports it so; ignore position which is incompatible with SNPedia info
      self.genotypes[tokens[0]] = Genotype(snp_name=tokens[0], genotype=tokens[3], chromosome=tokens[1], position=None, orientation='plus')
      if self.data.config.verbose and len(self.genotypes) % 100000 == 0 : 
        print 'INFO : loaded genome entry %d, %s' % (len(self.genotypes), tokens[0])
  
  def genotype(self, genotype_name) :
    try :
      return self.genotypes[genotype_name]
    except KeyError :
      return None
  
  def find_snps(self) :
    results = []
    for snp_name in self.data.snps() :
      print 'xxx Testing', snp_name
      if not snp_name in self.genotypes : continue
      genotype = self.genotypes[snp_name]
      snp = self.data.snp(snp_name)
      matches = snp.match(genotype)
      print 'xxx Test Results', ','.join(m.genotype for m in matches)
      if len(matches) == 0 : continue
      if len(matches) > 1 : 
        print 'WARNING : genotype %s of snp %s has multiple matches:' % (genotype.genotype, snp.name), ', '.join(g.genotype for g in matches)
        continue
      results.append(Match(genotype, snp, matches[0]))
    return results

  def find_genosets(self) :
    results = []
    for genoset_name in self.data.genosets() :
      print 'Genoset ' + genoset_name
      genoset = self.data.genoset(genoset_name)
      if genoset.match_genome(self) :
        print '==> passed', genoset_name
        results.append(genoset)
    return results


class Genotype :
  def __init__(self, snp_name, genotype, chromosome, position, orientation = 'plus') :
    self.snp_name = snp_name
    self.chromosome = chromosome
    self.position = position
    self.orientation = orientation
    self.genotype = genotype

  def alleles(self) :
    return Allele.split(self.genotype, self.orientation)


class Match :
  def __init__(self, genotype, snp, match) :
    self.genotype = genotype
    self.snp = snp
    self.snp_geno_match = match


class HtmlOutput :
  def __init__(self, config = None) :
    self.config = config if config != None else Config()
    
  def snp_page(self, snp_results, output, mag_cutoff = 2, mode = 'w') :
    print snp_results
    sorted_results = sorted(snp_results, key=lambda result: result.snp_geno_match.magnitude, reverse=True)
    outfile = open(output, mode)
    outfile.write(self.config.css_style)
    outfile.write('<h1> SNP Results </h1> <table>\n')
    outfile.write('<tr><td> SNP ID </td><td> Chr. </td><td> Gene </td><td> Individual Genotype </td><td> Matched Genotype </td><td> Significance </td><td> Description </td>\n')
    for result in sorted_results :
      if result.snp_geno_match.magnitude < mag_cutoff : break # list is sorted
      outfile.write('  <tr><td> %s </td><td> %s </td><td> %s </td> <td> %s </td>  <td> %s </td> ' 
                    % (result.snp.name, result.snp.chromosome, result.snp.gene, result.genotype.genotype, result.snp_geno_match.genotype))
      mag_style = ''
      if result.snp_geno_match.repute == 'good' : mag_style = '  bgcolor="#00FF00"'
      elif result.snp_geno_match.repute == 'bad' : mag_style = '  bgcolor="#FF0000"'
      outfile.write('<td%s> %s </td> <td> %s </td>' % (mag_style, result.snp_geno_match.magnitude, result.snp_geno_match.summary))
      outfile.write('<tr>\n')
    outfile.write('</table>\n')

  def genoset_page(self, genoset_results, output, mag_cutoff = 2, mode = 'w') :
    sorted_results = sorted(genoset_results, key=lambda gs: gs.magnitude, reverse=True)
    outfile = open(output, mode)
    outfile.write(self.config.css_style)
    outfile.write('<h1> Genoset Results </h1> <table>\n')
    outfile.write('  <tr><td> Genoset ID </td><td> Significance </td><td> Description </td>\n')
    for gs in sorted_results :
      if gs.magnitude < mag_cutoff : break # list is sorted
      outfile.write('  <tr><td> %s </td>' % gs.name)
      mag_style = ''
      #if results[4] == 'g' : mag_style = ' style="background-color:green'
      #elif results[4] == 'b' : mag_style = ' style="background-color:red'
      if gs.repute == 'good' : mag_style = '  bgcolor="#00FF00"'
      elif gs.repute == 'bad' : mag_style = '  bgcolor="#FF0000"'
      outfile.write('<td%s> %s </td> <td> %s </td>' % (mag_style, gs.magnitude, gs.summary))  # magnitude and summary
      outfile.write('<tr>\n')
    outfile.write('</table>\n')

  def snp_db_page(self, data, output, verbose = False, mode = 'w') :
    outfile = open(output, mode)
    outfile.write(css_style)
    #rs11155133 ['6', 'p', 'LOC102723724', '140848688', '38.1', '141', [('AA', ['0', 'common in complete genomics', 'g'])]]
    outfile.write('<h1> SNP Database </h1> <table>\n')
    if verbose :
      outfile.write('<tr><td> SNP </td> <td> Chromosome </td> <td> Orientation </td>  <td> Gene </td> <td> Location </td> <td> Genome Build </td>  <td> dbSNP Build </td> <td> Genotype </td> <td> Significance </td>  <td> Description </td></tr>\n')
    else:
      outfile.write('<tr><td> SNP </td> <td> Chromosome </td>  <td> Gene </td> <td> Genotype </td> <td> Significance </td>  <td> Description </td>\n')
    for snp_name in data.snps :
      snp = data.snps[snp_name]
      td = '<td rowspan="%d">' % len(snp.genotypes)
      if verbose :
        outfile.write('<tr>%s %s </td> %s %s </td>  %s %s </td> %s %s </td> %s %s </td>  %s %s </td> %s %s </td>' % 
                      (td, snp_name, td, snp.chromosome, td, snp.stabilized_orientation, td, snp.gene, td, snp.location, 
                       td, snp.genome_build, td, snp.dbSNP_build))
      else:
        outfile.write('<tr> %s %s </td> %s %s </td>  %s %s </td>' % 
                      (td, snp_name, td, snp.chromosome, td, snp.gene))
      for gt in snp.genotypes :
        mag_style = ''
        if gt.repute == 'g' : mag_style = ' bgcolor="#00FF00"'
        elif gt.repute == 'b' : mag_style = ' bgcolor="#FF0000"'
        outfile.write('<td> %s </td> <td%s> %s </td>  <td> %s </td></tr>\n' % (gt.name, mag_style, gt.magnitude, gt.summary))
    outfile.write('</table>\n')

  def genoset_db_page(self, data, genoset_db_file, output, mode = 'w') :
    outfile = open(output, mode)
    outfile.write(css_style)
    outfile.write('<h1> Genoset Database </h1> <table>\n')
    outfile.write('<tr><td> Genoset </td> <td> Significance </td> <td> Description </td> <td> Criteria </td></tr>\n')
    for genoset_name in data.genosets :
      genoset = data.genosets[genoset_name]
      mag_style = ''
      if genoset.repute == 'good' : mag_style = ' bgcolor="#00FF00"'
      elif genoset.repute == 'bad' : mag_style = ' bgcolor="#FF0000"'
      outfile.write('<td> %s </td> <td%s> %s </td>  <td> %s </td><td> %s </td></tr>\n' % 
                    (genoset_name, mag_style, genoset.magnitude, genoset.summary, genoset.criteria))
    outfile.write('</table>\n')

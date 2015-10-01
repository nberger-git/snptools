SNPTools : a set of tools to interpret genomes locally, without transmitting DNA data over the network.

The code relies on the SNPeida database (www.snpedia.com) to provide information on SNPs and "genosets", logical combinations on SNPs (and, or, etc.) to pinpoint certain conditions.

The processing has 2 main steps: first, the SNP and genoset data must be obtained from SNPedia. Seccond, this data is used to interpret the provided genome. Running the tool is therefore as simple as

cd snptools
./download.py  # download data from SNPedia
./interpret.py my_genome.txt output.html # interpret my_genome.txt and present results in output.html

However in practice the first step takes about 24h, since all ~76k SNPs must be downloaded one by one. The process can be interrupted however, and will resume where it left off when resumed.

The second step takes about 30 seconds, and produces a webpage summarizing the results. Only the 23andme genome format is supported for now.

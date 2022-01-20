#!/usr/bin/env python
from __future__ import division
import argparse
from itertools import islice

parser = argparse.ArgumentParser(description='Script to filter a sync file (Popoolation2) by a list of SNPs to keep (e.g. a snpdet file produced by poolfstat)')
parser.add_argument('-i', dest = 'sync', type = str, required=True,  help = 'input sync file')
parser.add_argument('-snps', dest = 'snps', type = str, required=True,  help = 'input snp file (Scaffold Position)')
parser.add_argument('-o', dest = 'out', type = str, required=True,  help = 'output filtered sync file')

args = parser.parse_args()

out_sync = open(args.out, 'w')

# list intersection function
def intersection(lst1, lst2): 
	lst3 = [value for value in lst1 if value in lst2] 
	return lst3 

## create a list of snps to keep 
snplist = []
with open(args.snps, 'rU') as f:
	for line in f:
		dat = line.split()
		snp = dat[0]+'\t'+dat[1]
		snplist.append(snp)
		
## turn the sync file into a dictionary with the SNP locations as the keys
## do this in chunks of 1 million snps to save memory
with open(args.sync) as f:
	while True:
		next_n_lines = list(islice(f, 1000000))
		if not next_n_lines:
			break		
		sync_dict = {}
		for line in next_n_lines:
			dat = line.split()
			snp = dat[0]+'\t'+dat[1]
			freqs = '\t'.join(dat[2:])
			sync_dict[snp] = freqs
		## find overlap	
		tokeep = intersection(snplist,sync_dict)

		## write out sync file for overlap
		for i in tokeep:
			outline = i + '\t' + sync_dict[i] + '\n'
			out_sync.write(outline)

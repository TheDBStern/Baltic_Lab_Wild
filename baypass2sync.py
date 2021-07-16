#!/usr/bin/env python
from __future__ import division
import glob, os
import sys
import argparse
import numpy

parser = argparse.ArgumentParser(description='Script to convert a file from multipopulation BayPass format including scaffold and position (e.g. scaffold position refcount1 altcount1 etc.) sync format')
parser.add_argument('-i', dest = 'input', type = str, required=True,  help = 'input file')
parser.add_argument('-out', dest = 'out', type = str, required=True,  help = 'sync output file')
parser.add_argument('-n', dest = 'numPops', type = int, required=True,  help = 'number of populations')


args = parser.parse_args()

sync_outfile = open(args.out,"w")

with open(args.input,'rU') as f:
	for line in f:
		chr=line.split(' ')[0]
		pos=line.split(' ')[1]
		pcount = 1
		genos = []
		for pop in range(0,args.numPops*2,2):
			c1 = line.split(' ')[pop+2] # get first allele count
			c2 = line.split(' ')[pop+3].strip('\n') # second allele count
			geno = "%s:%s:0:0:0:0"%(c1,c2)
			genos.append(geno)
			pcount +=1
		outgenos = '\t'.join(genos)
		sync_outfile.write(chr+'\t'+pos+'\t'+'A'+'\t'+outgenos+'\n')

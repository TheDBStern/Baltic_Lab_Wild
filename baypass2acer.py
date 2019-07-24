#!/usr/bin/env python
from __future__ import division
import glob, os
import sys
import argparse
import numpy

parser = argparse.ArgumentParser(description='Script to convert a file from multipopulation BayPass format (refcount1 altcount1 etc.) to frequencies of the alt allele and a coverage matrix')
parser.add_argument('-i', dest = 'input', type = str, required=True,  help = 'input file')
parser.add_argument('-freq', dest = 'freqout', type = str, required=True,  help = 'frequency output file')
parser.add_argument('-cov', dest = 'covout', type = str, required=True,  help = 'covarage output file')
parser.add_argument('-n', dest = 'numPops', type = int, required=True,  help = 'number of populations')
parser.add_argument('-fold', dest= 'fold', action ='store_true', default= False, help ='Calculate MAF for each SNP in each pop instead of ref allele frequency, default = False.')


args = parser.parse_args()

freq_outfile = open(args.freqout,"w")
cov_outfile = open(args.covout,"w")

with open(args.input,'rU') as f:
	for line in f:
		pcount = 1
		freqs = []
		covs = []
		for pop in range(0,args.numPops*2,2):
			c1 = int(line.split(' ')[pop]) # get first allele count for each pop
			c2 = int(line.split(' ')[pop+1]) # second allele count
			tot = c1 + c2
			if args.fold:
				freq = round(min(c1,c2)/tot,3)
			else:
				freq = round(c2/tot,3)
			freqs.append(str(freq))
			covs.append(str(tot))
			pcount +=1
		freq_outfile.write(' '.join(freqs)+'\n')
		cov_outfile.write(' '.join(covs)+'\n')
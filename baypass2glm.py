#!/usr/bin/env python

import glob, os
import argparse
import numpy

parser = argparse.ArgumentParser(description='Script to convert a file from BayPass format to one where allele counts are on two lines')
parser.add_argument('-i', dest = 'input', type = str, required=True,  help = 'input baypass format file')
parser.add_argument('-o', dest = 'output', type = str, required=True,  help = 'name of output')
parser.add_argument('-p', dest = 'pops', type = int, required=True,  help = 'number of pools')

args = parser.parse_args()

output = open(args.output, "w")

with open(args.input,'rU') as f:
	for line in f:
		A_counts = []
		a_counts = []
		for pop in range(0,args.pops*2,2):
			c1 = line.split(' ')[pop] # get first allele count for each pop
			c2 = line.split(' ')[pop+1] # second allele count
			A_counts.append(c1)
			a_counts.append(c2)
		output.write(' '.join(A_counts)+ '\n')
		output.write(' '.join(a_counts))
		

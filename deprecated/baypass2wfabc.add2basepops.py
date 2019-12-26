#!/usr/bin/env python

import glob, os
import argparse
import numpy

parser = argparse.ArgumentParser(description='Script to convert a file from BayPass format to WFABC format. Removes SNPs fixed in the base population. Assumes one shared base population which is the first population. Example: Base Rep1_T1 Rep1_T2 Rep2_T1 Rep2_T2...')
parser.add_argument('-i', dest = 'input', type = str, required=True,  help = 'input baypass format file')
parser.add_argument('-r', dest = 'reps', type = int, required=True,  help = 'number of replicates')
parser.add_argument('-nt', dest = 'timepoints', type = int, required=True,  help = 'number of timepoints sampled beyond the base')
parser.add_argument('-gen', dest = 'gen', nargs='+', type = str, required=True,  help = 'generation number for each timepoint')

args = parser.parse_args()

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
flen = file_len(args.input)

#create file and headers for each rep
for i in range(1,args.reps+1):
	output = open('rep'+str(i)+'.snps','w')
	output.write(str(flen)+' '+str(args.timepoints+1)+'\n')
	output.write('0'+','+','.join(args.gen)+'\n')

#create excluded snp files
excluded = open('excludedSNPs.txt','w')

with open(args.input,'rU') as f:
	lcount = 1
	for line in f:
		pcount = 1
		base_c1a = int(line.split(' ')[0])
		base_c1b = int(line.split(' ')[2])
		base_c2a = int(line.split(' ')[1])
		base_c2b = int(line.split(' ')[3])
		base_tot = base_c1a + base_c1b + base_c2a + base_c2b
		base_A = base_c1a + base_c1b
		base_a = base_c2a + base_c2b
		if base_A == 0 or base_a == 0:
			excluded.write(str(lcount)+'\n')
		else:
			for rep in range(4,args.reps*2*args.timepoints+4,args.timepoints*2):
				repdat = line.split(' ')[rep:rep+args.timepoints*2]
				totals = []
				a_count = []
				for i in range(0,args.timepoints*2,2):
					c1 = repdat[i] # get first allele count
					c2 = repdat[i+1] # second allele count
					tot = str(int(c1) + int(c2))
					totals.append(tot)
					a_count.append(c2)
				output = open('rep'+str(pcount)+'.snps','a')
				output.write(str(base_tot)+','+','.join(totals)+'\n'+str(base_a)+','+','.join(a_count)+'\n')
				pcount +=1
		lcount +=1

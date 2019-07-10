#!/usr/bin/env python

import glob, os
import argparse
import numpy

parser = argparse.ArgumentParser(description='Script to convert a file from BayPass format to WFABC format. Assumes one shared base population which is the first population')
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
	output = open('rep'+str(i),'w')
	output.write(str(flen)+' '+str(args.timepoints+1)+'\n')
	output.write('0'+','+','.join(args.gen)+'\n')

with open(args.input,'rU') as f:
	for line in f:
		pcount = 1
		for rep in range(2,args.reps*2*args.timepoints+2,args.timepoints*2):
			base_c1 = int(line.split(' ')[0])
			base_c2 = int(line.split(' ')[1])
			base_tot = base_c1 + base_c2
			repdat = line.split(' ')[rep:rep+args.timepoints*2]
			totals = []
			A_count = []
			for i in range(0,args.timepoints*2,2):
				c1 = repdat[i] # get first allele count
				c2 = repdat[i+1] # second allele count
				tot = str(int(c1) + int(c2))
				totals.append(tot)
				A_count.append(c2)
			output = open('rep'+str(pcount),'a')
			output.write(str(base_tot)+','+','.join(totals)+'\n'+str(base_c2)+','+','.join(A_count)+'\n')
			pcount +=1


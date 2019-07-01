#!/usr/bin/env python

import argparse


parser = argparse.ArgumentParser(description='Script to filter a sync file (Popoolation2) by a list of SNPs to keep (e.g. a snpdet file produced by poolfstat)')
parser.add_argument('-i', dest = 'sync', type = str, required=True,  help = 'input sync file')
parser.add_argument('-snps', dest = 'snps', type = str, required=True,  help = 'input snp file (Scaffold Position)')
parser.add_argument('-o', dest = 'out', type = str, required=True,  help = 'output filtered sync file')

args = parser.parse_args()

out_sync = open(args.out, 'w')

## create a list of snps to keep 
snplist = []
with open(args.snps, 'rU') as f:
	for line in f:
		dat = line.split()
		snp = dat[0]+' '+dat[1]
		snplist.append(snp)

## go through sync file line by line and keep only snps in snp list
with open(args.sync, 'rU') as f:
	for line in f:
		dat = line.split()
		snp = dat[0]+' '+dat[1]
		if snp in snplist:
			out_sync.write(line)
		else:
			pass
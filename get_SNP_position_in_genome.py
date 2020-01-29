#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description='Script to convert pseudoreference SNP position to approximate position in genome based on best blast hit')
parser.add_argument('-i', dest = 'snps', type = str, required=True,  help = 'input snp file, format `transcript position`')
parser.add_argument('-b', dest = 'blast', type = str, required=True,  help = 'input blast output file, outfmt 6')
parser.add_argument('-o', dest = 'out', type = str, required=True,  help = 'output converted SNP file')
args = parser.parse_args()

out_snps = open(args.out, 'w')
blastfile = open(args.blast, 'rU')

with open(args.snps,'rU') as f:
	for l in f:
		contig = l.split()[0]
		position = int(l.split()[1])
		with open(args.blast, 'rU') as b:
			for line in b:
				# the first entry is always the best hit based on bitscore
				if line.split('\t')[0] == contig:
					scaffold = line.split('\t')[1]
					genomePosition1 = int(line.split('\t')[8])
					genomePosition2 = int(line.split('\t')[9])
					transPosition1 = int(line.split('\t')[6])
					transPosition2 = int(line.split('\t')[7])
					# an approximate position for the SNP could be obtained by finding the approximate start position and add the SNP position
					if transPosition1 < transPosition2:
						start = transPosition1
					else:
						start = transPosition2				
					if genomePosition1 < genomePosition2:
						snpPosition =  genomePosition1 - start + position
					else:
						snpPosition =  genomePosition2 - start + position
					out_snps.write(scaffold+'\t'+str(snpPosition)+'\n')
					b.close()
					break
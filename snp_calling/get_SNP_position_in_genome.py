#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description='Script to convert SNP positions called in one reference genome to approximate position in another genome based on blast results')
parser.add_argument('-i', dest = 'snps', type = str, required=True,  help = 'input snp file, format `transcript position`')
parser.add_argument('-b', dest = 'blast', type = str, required=True,  help = 'input blast output file, outfmt 6')
parser.add_argument('-o', dest = 'out', type = str, required=True,  help = 'output converted SNP file')
args = parser.parse_args()

out_snps = open(args.out, 'w')

### create dictionaries
scaffold_dict = {}
genomePosition1_dict = {}
genomePosition2_dict = {}
transPosition1_dict = {}
transPosition2_dict = {}

print('Creating blast output dictionary')
current_gene = ''
with open(args.blast, 'rU') as b:
	for line in b:
		gene = line.split('\t')[0]
		# the first entry is always the best hit based on bitscore
		if len(current_gene) == 0: # deal with first line
			current_gene = gene
			scaffold_dict[gene] = line.split('\t')[1]
			genomePosition1_dict[gene] = int(line.split('\t')[8])
			genomePosition2_dict[gene] = int(line.split('\t')[9])
			transPosition1_dict[gene] = int(line.split('\t')[6])
			transPosition2_dict[gene] = int(line.split('\t')[7])
		elif current_gene == gene:
			pass
		else:
			current_gene = gene
			scaffold_dict[gene] = line.split('\t')[1]
			genomePosition1_dict[gene] = int(line.split('\t')[8])
			genomePosition2_dict[gene] = int(line.split('\t')[9])
			transPosition1_dict[gene] = int(line.split('\t')[6])
			transPosition2_dict[gene] = int(line.split('\t')[7])

print('Blast dictionary finished')

with open(args.snps,'rU') as f:
	for line in f:
		contig = line.split()[0]
		position = int(line.split()[1])
		if contig in scaffold_dict.keys():
			scaffold = scaffold_dict[contig]
			genomePosition1 = genomePosition1_dict[contig]
			genomePosition2 = genomePosition2_dict[contig]
			transPosition1 = transPosition1_dict[contig]
			transPosition2 = transPosition2_dict[contig]
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
		else:
			out_snps.write('NA'+'\t'+'NA'+'\n')

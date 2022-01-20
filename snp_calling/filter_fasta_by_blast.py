#!/usr/bin/env python
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Script to filter a multifasta file based on whether sequences had a significant blast hit to some database')
parser.add_argument('-i', dest = 'fasta', type = str, required=True,  help = 'input fasta file')
parser.add_argument('-b', dest = 'blast', type = str, required=True,  help = 'input blast output file, outfmt 6')
parser.add_argument('-o', dest = 'out', type = str, required=True,  help = 'output filtered fasta file')
parser.add_argument('-e', dest = 'evalue', type = float, default=1e-3,  help = 'e-value cutoff, default 1e-3')
args = parser.parse_args()


in_fasta = SeqIO.index(args.fasta, 'fasta')
out_fasta = open(args.out, 'w')

current_gene = ''
evals = []
out_genes = []
with open(args.blast,'rU') as f:
	for line in f:
		gene = line.split('\t')[0]
		if len(current_gene) == 0: # deal with first line
			current_gene = gene
			evals.append(float(line.split('\t')[10]))
		elif current_gene == gene:
			evals.append(float(line.split('\t')[10]))
		else:
			if any(x <= args.evalue for x in evals):
				if gene not in out_genes:
					out_genes.append(gene)
			current_gene = gene
			evals = []
			evals.append(float(line.split('\t')[10]))

print(str(len(out_genes))+' genes passed evalue filter')

for gene in out_genes:
	SeqIO.write(in_fasta[gene], out_fasta, 'fasta')
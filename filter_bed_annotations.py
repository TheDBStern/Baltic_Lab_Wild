#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description='Script to filter gtf based on whether genes have GO terms assigned')
parser.add_argument('-i', dest = 'input', type = str, required=True,  help = 'input bed')
parser.add_argument('-go', dest = 'goterms', type = str, required=True,  help = 'go term file. Format: GoNumber\tDescription\tGeneIDs')
parser.add_argument('-o', dest = 'output', type = str, required=True,  help = 'output filtered bed')
args = parser.parse_args()


genes = []

with open(args.goterms, 'rU') as f:
	for line in f:
		geneids = line.split('\t')[2]
		for gene in geneids.split(' '):
			if gene not in genes:
				genes.append(gene)

print('%s genes have GO Terms'%(len(genes)))
				
outfile = open(args.output, 'w')

with open(args.input, 'rU') as f:
	for line in f:
		gene_id = line.split('\t')[3]
		if gene_id in genes:
			outfile.write(line)
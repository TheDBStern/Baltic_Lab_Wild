#!/usr/bin/env python
import argparse
from BCBio import GFF

parser = argparse.ArgumentParser(description='Script to filter gtf based on whether genes have GO terms assigned')
parser.add_argument('-i', dest = 'input', type = str, required=True,  help = 'input gtf')
parser.add_argument('-go', dest = 'goterms', type = str, required=True,  help = 'go term file. Format: GoNumber\tDescription\tGeneIDs')
parser.add_argument('-o', dest = 'output', type = str, required=True,  help = 'output filtered gtf')
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
		dat = line.split('\t')[8]
		for item in dat.split(';'):
			info = item.lstrip(' ')
			#print(info)
			if info.split(' ')[0] == "gene_id":
				gene_id = info.split(' ')[1].strip('"')
				print(gene_id)
				if gene_id in genes:
					outfile.write(line)
#!/usr/bin/env python

import sys, os
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Create a .gtf file from a fasta file, considering each scaffold a gene')
parser.add_argument('-i', dest= 'fasta', type = str, required=True, help ='Input fasta file')
parser.add_argument('-o', dest= 'gtf', type = str, required=True, help ='Output gtf file')
args = parser.parse_args()

output = open(args.gtf, 'w')

fasta = SeqIO.parse(args.fasta, 'fasta')

for record in fasta:
	outline = '%s\tpseudoref\texon\t1\t%s\t.\t+\t.\tgene_id "%s";\n'%(record.id,len(record.seq),record.id)
	output.write(outline)

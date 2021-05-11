#!/usr/bin/env python

import sys, os
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='For a set of left/R1 reads, fetch corresponding right/R2 read pairs')
parser.add_argument('-m', dest= 'mapped', type = str, required=True, help ='Mapped left reads file in fastq format for which you want to get corresponding right reads')
parser.add_argument('-r', dest= 'right', type = str, required=True, help ='Full right reads file in fastq format')
parser.add_argument('-o', dest= 'out', type = str, required=True, help ='Output corresponding right reads')
args = parser.parse_args()

print 'Indexing reads...'
right_idx = SeqIO.index(args.right, 'fastq')

output = open(args.out, 'w')
left = SeqIO.parse(args.mapped, 'fastq')
for seq in left:
	if seq.id.endswith('/1'):
		right_id = seq.id.replace('/1','/2')
		SeqIO.write(right_idx[right_id], output, 'fastq')
	else:
		right_id = seq.id.replace(' 1:',' 2:')
		SeqIO.write(right_idx[right_id], output, 'fastq')

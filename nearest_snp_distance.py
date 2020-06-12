#! /user/local/bin/python

import argparse

parser = argparse.ArgumentParser(description='For every SNP, calculates distance to nearest neighboring SNP.')
parser.add_argument('-i', dest = 'snps', type = str, required=True,  help = 'input snp file. Format `scaffold position`. Must be sorted by scaffold and position')
parser.add_argument('-o', dest = 'out', type = str, required=True,  help = 'output file')
args = parser.parse_args()

output = open(args.out,'w')

back_dist = ''
for_dist = ''
last_scaf = ''
last_pos = ''
with open(args.snps,'r') as f:	
	for line in f:
		scaf = line.split()[0]
		position = int(line.split()[1])
		if len(last_scaf) == 0: #deal with first line
			last_scaf = scaf
			last_pos = position
		elif scaf == last_scaf: #scaffold is same as previous, calculate distance and store. Check which dist is shorter
			for_dist = position - last_pos
			if isinstance(back_dist,int) and for_dist > back_dist:
				output.write(last_scaf+'\t'+str(last_pos)+'\t'+str(back_dist)+'\n')
			elif isinstance(back_dist,int) and for_dist < back_dist:
				output.write(last_scaf+'\t'+str(last_pos)+'\t'+str(for_dist)+'\n')
			elif not isinstance(back_dist,int): # first snp in scaffold
				output.write(last_scaf+'\t'+str(last_pos)+'\t'+str(for_dist)+'\n')
			back_dist = for_dist
			last_scaf = scaf
			last_pos = position
		elif scaf != last_scaf: #scaffold is not the same as previous SNP 
			if not isinstance(back_dist,int) and not isinstance(for_dist,int): #only one snp on scaffold
				pass
			else:
				output.write(last_scaf+'\t'+str(last_pos)+'\t'+str(back_dist)+'\n')
				last_scaf = scaf
				last_pos = position
				back_dist = ''
				for_dist = ''
'''
to do: how to deal with last line in the file
'''
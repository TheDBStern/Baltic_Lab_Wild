#!/usr/bin/env python

import sys

toflip = open(sys.argv[1],'rU')
output = open(sys.argv[2],'w')

godict = {}


for line in toflip:
	trans = line.split('\t')[0]
	goterms = line.split('\t')[1].split(',')
	for term in goterms:
		if term in godict:
			tmplist = godict[term]
			tmplist.append(trans)
			godict[term] = tmplist
		else:
			godict[term] = [trans]

#print(godict)
for key in godict.keys():
	goterm=key
	print(goterm)
	with open('/Users/dbstern/Desktop/Genomics_Programs/go-basic.obo','rU') as f:
		for lines in f:
			if "id: "+goterm == lines.strip('\n'):
				info = next(f)
				description = info.split(':')[1].strip('\n')
				print(description)
				output.write(goterm+'\t'+description+'\t'+' '.join(godict[goterm])+'\n')	
		f.close()

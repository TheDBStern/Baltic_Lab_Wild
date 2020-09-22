#!/bin/bash
input="sig_snps.cmh05.lrt01.transcripts.txt"
while IFS= read -r line
do
	gunzip -c /Users/dbstern/Desktop/Baltic_sea_project/pseudoref/reference/trinotate_annotation_report.xls.gz | grep "$line"
	#echo "$line"
done < "$input"
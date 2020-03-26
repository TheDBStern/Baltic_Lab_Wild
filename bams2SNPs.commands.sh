## Commands to call SNPs from multiple bam files, creating VCF, genobaypass and sync files

samtools mpileup -q 20 -Q 0 -d 8000 -R -A –B -f Iteration_2.trinity.Trinity.cdhit95.filtered.fasta -b bamfiles.txt -o lab.mpileup

## Generate .vcf and genobaypass files along with pairwise Fst using poolfstat
varscan mpileup2cns lab.mpileup --min-coverage 20 --min-avg-qual 20 --min-reads2 4 --min-var-freq 0.001 --output-vcf 1 --variants > lab.vcf

## max depth per pop is based on overall 99% - 99.9% read count distribution (calculate_coverage_distribution_sync.py)
Rscript vcf2genobaypass.R

## Generate sync file using Popoolation

# Create the sync file, filtering bases with quality less than 20
java -ea -Xmx50g -jar <path_to_popoolation2>/mpileup2sync.jar --input lab.mpileup --output lab.sync --min-qual 20 --threads 8

### Filter the sync file to match the SNPs called with poolfstat
python filter_sync_by_snplist.py -i lab.sync -snps lab.snpdet -o lab.filtered.sync
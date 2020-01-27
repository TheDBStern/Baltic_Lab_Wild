## Commands to call SNPs from multiple bam files, creating VCF, genobaypass and sync files


########
## Lab samples
########

##bamlist.txt (SNP frequencies will be in this order)
#BSE-0.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-1-T1.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-1-T2.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-2-T1.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-2-T2.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-3-T1.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-3-T2.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-4-T1.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-4-T2.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-5-T1.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-5-T2.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-6-T1.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-6-T2.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-8-T1.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-8-T2.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-9-T1.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-9-T2.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-11-T1.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-11-T2.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-12-T1.ngm.mapq3.sorted.rg.dedup.real.bam
#BSE-12-T2.ngm.mapq3.sorted.rg.dedup.real.bam


samtools mpileup -q 3 -Q 0 -d 8000 -R -A –B -f Eaff_11172013.genome.masked.fa -b bamfiles.txt -o lab.base_T1_T2.mpileup

## Generate .vcf and genobaypass files along with pairwise Fst using poolfstat
varscan mpileup2cns lab.base_T1_T2_10.mpileup --p-value 0.99 --min-coverage 50 --min-avg-qual 20 --min-reads2 4 --min-var-freq 0.001 --output-vcf 1 --variants > lab.base_T1_T2.vcf

## max depth per pop is based on overall 99% - 99.9% read count distribution (calculate_coverage_distribution_sync.py)
Rscript vcf2genobaypass.R

## Generate sync file using Popoolation

# Create the sync file, filtering bases with quality less than 20
java -ea -Xmx50g -jar <path_to_popoolation2>/mpileup2sync.jar --input lab.base_T1_T2.mpileup --output lab.base_T1_T2.sync --min-qual 20 --threads 8


### Filter the sync file to match the SNPs called with poolfstat
python filter_sync_by_snplist.py -i lab.base_T1_T2.sync -snps lab.base_T1_T2.snpdet -o lab.base_T1_T2.filtered.sync
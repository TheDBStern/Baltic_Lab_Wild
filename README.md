# Introduction
This repository contains analysis scripts and data associated with our manuscript

> Stern DB, Anderson NW, Diaz JA, and CE Lee. Parallel polygenic adaptation promoted by epistasis in laboratory and wild populations of a Baltic Sea Copepod

## Usage
Scripts are organized by either [**snp_calling**](./snp_calling) or [**selection_analyses**](./selection_analyses)  
Command-line options for python scripts can be found, e.g.,`baypass2freqs_cov.py -h`

**snp_calling**, i.e., reference assembly, SNP calling, SNP data processing
- [assemble_poolseq.commands.txt](./snp_calling/assemble_poolseq.commands.txt) Commands used to generate the 'pseudoreference' genome by 'tiling' Pool-seq data onto the transcriptome in an iterative mapping and assembly approach
- [baypass2freqs_cov.py](./snp_calling/baypass2freqs_cov.py) Converts a file from multipopulation BayPass format (refcount1 altcount1 etc.) to frequencies of the alt allele and a coverage matrix
- [bams2SNPs.commands.sh](./snp_calling/bams2SNPs.commands.sh) Commands used to call SNPs and generate allele count files
- [calculate_coverage_distribution_sync.py](./snp_calling/calculate_coverage_distribution_sync.py) Calculates the top X percentage of coverage across all pools from a sync file
- [filter_fasta_by_blast.py](./snp_calling/filter_fasta_by_blast.py) Filters a multifasta file based on whether sequences had a significant blast hit to some sequence database or genome
- [filter_sync_by_snplist.py](./snp_calling/filter_sync_by_snplist.py) Filters a sync file (Popoolation2) by a list of SNPs to keep (e.g. a snpdet file produced by poolfstat)
- [get_mates.py](./snp_calling/get_mates.py) For a set of left/R1 reads, fetch corresponding right/R2 read pairs
- [get_SNP_position_in_genome.py](./snp_calling/get_SNP_position_in_genome.py) Convert SNP positions called in one reference genome to approximate position in another genome based on blast results
- [vcf2genobaypass.R](./snp_calling/vcf2genobaypass.R) R commands to generate the read count file from the VarScan VCF using *poolfstat*

**selection_analyses**, i.e., CMH, Chi-square, & LMM tests, calculating Jaccard index
- [ACER_code.R](./selection_analyses/ACER_code.R) R commands used to run the Chi-square and CMH tests on SNPs
- [determine_AFC_cutoff.R](./selection_analyses/determine_AFC_cutoff.R) R commands to simulate neutral allele frequency change to determine a cutoff to call an allele an under selection in a given line
- [parallelism_functions.R](./selection_analyses/parallelism_functions.R) R functions to calculate the Jaccard index and RFS for the empirical data
- [prep_lmm.R](./selection_analyses/prep_lmm.R) R code specific to this study for generating the input file to run the lmm analysis of SNP frequency trajectories. Uses the files in the [data](./data) directory
    * 'prep_lmm.rawAFC.R' - same as above but does not transform the allele frequencies
    * 'prep_lmm.rawFreqs.R' - same as above but uses raw allele frequencies rather than divergence from the ancestor
- [run_lmm.R](./selection_analyses/run_lmm.R) R script to run the linear mixed model with lme4 on every called SNP. Uses the output from [prep_lmm.R](./prep_lmm.R)

## Software required to run these scripts
- [BLAST 2.7.1+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [BWA-MEM v0.7.17](http://bio-bwa.sourceforge.net/bwa.shtml)
- [CD-HIT v4.7](http://weizhongli-lab.org/cd-hit/)
- [PoPoolation2](https://sourceforge.net/p/popoolation2/wiki/Main/)
- [SAMBLASTER v0.1.26](https://github.com/GregoryFaust/samblaster)
- [Samtools v1.3.1](http://www.htslib.org/)
- [Trinity v2.6.6](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
- [VarScan v2.4.3](http://varscan.sourceforge.net/)

## Python packages
Python version 3.8.2
- [BioPython v1.78](https://biopython.org/)
- [joblib v.1.0.1](https://joblib.readthedocs.io/en/latest/)
- [numpy v1.15.2](https://numpy.org/)

## R packages
R version 4.0.4
- [data.table v.1.14.0](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html)
- [dplyr v1.0.5](https://dplyr.tidyverse.org/)
- [lme4 v1.1.21](https://cran.r-project.org/web/packages/lme4/lme4.pdf)
- [poolfstat v1.1.1](https://cran.r-project.org/web/packages/poolfstat/poolfstat.pdf)
- [qvalue v2.14.1](https://github.com/StoreyLab/qvalue)
- [tcR v2.3.2](https://cran.r-project.org/web/packages/tcR/index.html)

## Other software used in the manuscript
- [ACER v1.0.2](https://github.com/MartaPelizzola/ACER)
- [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/)
- [BEDOPS v2.4.39](https://bedops.readthedocs.io/en/latest/)
- [Bowtie v2.3.5](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [DAVID v6.8](https://david.ncifcrf.gov/)
- [haplovalidate v0.1.4](https://github.com/kathrinannaotte/haplovalidate)
- [HMMER v3.2.1](http://hmmer.org/)
- [RSEM v1.3.1](https://deweylab.github.io/RSEM/)
- [PolygenicAdaptationCode](https://github.com/jjberg2/PolygenicAdaptationCode)
- [Transdecoder v5.5](https://github.com/TransDecoder/TransDecoder/wiki)
- [Trimmomatic v0.39](http://www.usadellab.org/cms/?page=trimmomatic)
- [TreeMix v1.13](https://bitbucket.org/nygcresearch/treemix/wiki/Home)

## Data
SNPs and allele counts derived from the Pool-seq data are available in the [data](./data) directory. Please see the README file within for information.

Please contact the authors for questions or issues.

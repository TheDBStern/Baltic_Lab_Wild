## Data description

### lab.snps
**lab.snps_freq** and **lab.snps_cov** contain SNP frequencies and read coverage for every SNP in every laboratory line,
generated using R::*poolfstat* and *baypass2freqs_cov.py*  
**lab.snpdet** contains the reference scaffold and position for every SNP
in the SNP frequency and coverage files

### hap_blocks
**hap_blocks.freq** and **hap_blocks.cov** contain the frequencies and coverages (based on the median SNP frequency in each haplotype blocks)
for all 121 haplotype blocks in all  
**hap_blocks.res.RDS** is an rds file (read into R using readRDS(hap_blocks.res.RDS)) with information for each haplotype block (start and stop position, starting frequency, selection coef., etc.)

The data columns are in the following order:

| Sample ID | Description |
| --- | --- |
| BSE-0 | Starting population, sample 1 |
| BSE-0B | Starting population, sample 2 |
| BS3C | Control line 3C, generation 20 |
| BS4C | Control line 4C, generation 20 |
| BSE-1-T1 | Control line 1, generation 6 |
| BSE-1-T2 | Control line 1, generation 10 |
| BSE-1-T2B | Control line 1, generation 10, technical replicate |
| BSE-2-T1 | Control line 2, generation 6 |
| BSE-2-T2 | Control line 22, generation 10 |
| BSE-2-T2B | Control line 2, generation 10, technical replicate |
| BSE-3-T1 | Treatment line 3, generation 6 |
| BSE-3-T2 | Treatment line 3, generation 10 |
| BSE-4-T1 | Treatment line 4, generation 6 |
| BSE-4-T2 | Treatment line 4, generation 10 |
| BSE-5-T1 | Treatment line 5, generation 6 |
| BSE-5-T2 | Treatment line 5, generation 10 |
| BSE-6-T1 | Treatment line 6, generation 6 |
| BSE-6-T2 | Treatment line 6, generation 10 |
| BSE-7-T1 | Treatment line 7, generation 6 |
| BSE-8-T1 | Treatment line 8, generation 6 |
| BSE-8-T2 | Treatment line 8, generation 10 |
| BSE-9-T1 | Treatment line 9, generation 6 |
| BSE-9-T2 | Treatment line 9, generation 10 |
| BSE-10-T1 | Treatment line 10, generation 6 |
| BSE-11-T1 | Treatment line 11, generation 6 |
| BSE-11-T2 | Treatment line 11, generation 6 |
| BSE-12-T1 | Treatment line 12, generation 6 |
| BSE-12-T2 | Treatment line 12, generation 6 |

### wild.snps and wild.snpdet
**wild.snps_freq** and **wild.snps_cov** contains SNP frequencies and coverages for every SNP in every wild population,
generated using R::*poolfstat* and *baypass2freqs_cov.py*  
**wild.snpdet** contains the reference scaffold and position for every SNP
in the SNP frequency and coverage files

The data columns are in the following order:

| Sample ID | Description |
| --- | --- |
| BB1E | Bothnian Bay 1 |
| BB2E | Bothnian Bay 2  |
| GBE | Gulf of Bothnia  |
| HF1E | Helsinki Field Station 1 |
| HF2E | Helsinki Field Station 2 |
| KIE | Kiel, Germany |
| RG1E | Gulf of Riga 1 |
| RG2E | Gulf of Riga 2|
| STE | Stockholm, Sweden |

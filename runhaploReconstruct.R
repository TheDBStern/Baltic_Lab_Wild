## Commands to run haploReconstruct on selection experiment with T0, T1, and T2

##import haploReconstruct library
library(haploReconstruct)

## Create variable that formats sync file into haploReconstruct format with an added header
#The file contains samples for T0, T1, and T2, each for eight replicate simulations
data=sync_to_frequencies(file="lab.base_T1_T2_noctrl.haprecon.sync", base.pops=rep(c(TRUE, rep(FALSE,2)),times=8), header=FALSE)

## Filter replicated time series data for informative SNPs
#default values used for ***put in list of parameters***, winsize=100KB
data_filtered=initialize_SNP_time_series(chr=data$chr, pos=data$pos, base.freq=data$basePops, lib.freqs=data[,7:ncol(data)], pop.ident=rep(1:8, each=3), pop.generation=rep(c(0,6,10), times=8), use.libs=rep(TRUE,24), winsize=100000)

## Reconstruct haplotype-blocks
#default values used for ***put in list of parameters***
#create list of unique scaffold values, then reconstruct haplotype blocks on each unique scaffold
Scaff_in_Filt <- data_filtered@col.info
Scaffolds <- unique(data_filtered$chr)
Scaffolds <- as.character(Scaffolds)
data_reconst <- lapply(Scaffolds, function(x) reconstruct_hb(data_filtered, chrom=x))

## Save data_reconst object as .RDS file on computer
saveRDS(data_reconst, "hapreconst_winsize100KB.RDS")

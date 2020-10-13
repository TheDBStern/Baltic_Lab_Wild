### Generate mimicree2 haplotype file for selected SNPs assuming each SNP is on a different chromosome (i.e. unlinked)

library(dplyr)
# get starting allele frequencies
sig <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.new.cmh05.lrt05.RDS')
afs <- sig$T0_AF_rising

#set effective pop size
ne <- 2000

genotypes <- c()

for (i in 1:length(afs)){
  af <- afs[i]
  genos <- rbinom(2*ne,1,af)
  genos[which(genos==0)] <- 'T'
  genos[which(genos==1)] <- 'A'
  genos_h1 <- genos[1:ne]
  genos_h2 <- genos[(ne+1):(2*ne)]
  hap <- paste(genos_h1,genos_h2,sep='',collapse=' ')
  genotypes <- c(genotypes,hap)
}

hap_out <- data.frame(1:length(afs),rep(10000,length(afs)),rep('A',length(afs)),rep('A/T',length(afs)),genotypes)
write.table(hap_out,"mimicree2_haplotypes.mimhap",sep='\t',quote=F,row.names=F,col.names = F)

## Generate effect size file based on selection coefficients
effs <- round(abs(sig$selCoef),4)

eff_out <- data.frame(1:length(effs),rep(10000,length(effs)),rep('A/T',length(effs)),effs,effs/2)
write.table(eff_out,"effect_sizes.txt",sep='\t',quote=F,row.names=F,col.names = F)

### figure out fitness function
library(ggplot2)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

fun.dir <- function(x) 0.5 + 500/(1+exp((-1*x)+15))
p + stat_function(fun = fun.dir) + xlim(0,15) + ylim(0,100)

fun.dim <- function(x) 0.0 + 1.2*(1-1/exp(0.02*x))
p + stat_function(fun = fun.dim) + xlim(-5,20) + ylim(0,1.2)

fun.exp <- function(x) 0.5 + exp(x)
p + stat_function(fun = fun.exp) + xlim(-5,10) + ylim(0,20)

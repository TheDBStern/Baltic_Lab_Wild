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

fun.dir <- function(x) 0.5 + 500/(1+exp((-1*x)+5))
p + stat_function(fun = fun.dir) + xlim(0,15) + ylim(0,500)

fun.dim <- function(x) 0.0 + 1.2*(1-1/exp(0.02*x))
p + stat_function(fun = fun.dim) + xlim(-5,20) + ylim(0,1.2)

fun.exp <- function(x) 0.5 + exp(x)
p + stat_function(fun = fun.exp) + xlim(-5,10) + ylim(0,20)

##prep interpolated fitness function as exponential
#average phenotype of starting population is ~-7.68
p <- seq(-10,20,0.1)
#f <- exp(1.5*p/7.68) +1
f <- (p+8.68)^3
#f <- exp((p+8.68/1000))
plot(p,f)

df <- data.frame(rep(1,length(p)), p, f)

write.table(df,'ff_interp.pos_epi.txt',quote=F,row.names=F,col.names=F,sep='\t')

## recombination simulations for 10Mb genome
library(dplyr)
# get starting allele frequencies
all <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.new.RDS')

afs <- sample(all$T0_AF_rising[which(all$T0_AF_rising != 0 & all$T0_AF_rising != 1)],10000)

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

positions <- sample(1:9999999, 10000)
hap_out <- data.frame("chr"=rep("Ch1",length(afs)),"pos"=positions,"ref"=rep('A',length(afs)),"alt"=rep('A/T',length(afs)),"geno"=genotypes)

g1 <- c('A',rep('T',ne-1))
g2 <- rep('T',ne)
geno <- paste(g1,g2,sep='',collapse=' ')
selected_snp <- data.frame("chr"="Ch1","pos"="5000000","ref"='A',"alt"='A/T',"geno"=geno)

hap_out <- rbind(hap_out,selected_snp)
hap_out <- hap_out %>%
          arrange(chr,pos)
write.table(hap_out,"haplotypes.10Mb.sel5M.ne2000.mimhap",sep='\t',quote=F,row.names=F,col.names = F)


## simulations for whole genome (4 chromsomes, 400 Mb), all selected SNPs, with recombination
all <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.new.RDS')
afs <- all$T0_AF[which(all$T0_AF_rising != 0 & all$T0_AF_rising != 1)]

sig <- readRDS('lab.res.new.cmh05.lrt05.RDS')
sig_afs <- sig$T0_AF[which(sig$T0_AF != 0 & sig$T0_AF != 1)]
#set effective pop size
ne <- 2000

get_geno <- function(af,ne){
  genos <- rbinom(2*ne,1,af)
  genos[which(genos==0)] <- 'T'
  genos[which(genos==1)] <- 'A'
  genos_h1 <- genos[1:ne]
  genos_h2 <- genos[(ne+1):(2*ne)]
  hap <- paste(genos_h1,genos_h2,sep='',collapse=' ')
  return(hap)
}

Ch_length <- round((length(afs)-length(sig_afs))/4,0)

Ch1_positions <- sample(1:99999999, Ch_length)
Ch1_sig_positions <- sample(Ch1_positions,round(length(sig_afs)/4))
Ch1_neut_positions <- Ch1_positions[which(!(Ch1_positions %in% Ch1_sig_positions))]
Ch1_geno <- sapply(sample(afs,length(Ch1_neut_positions)), get_geno, ne=2000)
Ch1_sig_geno <- sapply(sample(sig_afs,length(Ch1_sig_positions)), get_geno, ne=2000)

Ch2_positions <- sample(1:99999999, Ch_length)
Ch2_sig_positions <- sample(Ch2_positions,round(length(sig_afs)/4))
Ch2_neut_positions <- Ch2_positions[which(!(Ch2_positions %in% Ch2_sig_positions))]
Ch2_geno <- sapply(sample(afs,length(Ch2_neut_positions)), get_geno, ne=2000)
Ch2_sig_geno <- sapply(sample(sig_afs,length(Ch2_sig_positions)), get_geno, ne=2000)

Ch3_positions <- sample(1:99999999, Ch_length)
Ch3_sig_positions <- sample(Ch3_positions,round(length(sig_afs)/4))
Ch3_neut_positions <- Ch3_positions[which(!(Ch3_positions %in% Ch3_sig_positions))]
Ch3_geno <- sapply(sample(afs,length(Ch3_neut_positions)), get_geno, ne=2000)
Ch3_sig_geno <- sapply(sample(sig_afs,length(Ch3_sig_positions)), get_geno, ne=2000)

Ch4_positions <- sample(1:99999999, Ch_length)
Ch4_sig_positions <- sample(Ch4_positions,round(length(sig_afs)/4))
Ch4_neut_positions <- Ch4_positions[which(!(Ch4_positions %in% Ch4_sig_positions))]
Ch4_geno <- sapply(sample(afs,length(Ch4_neut_positions)), get_geno, ne=2000)
Ch4_sig_geno <- sapply(sample(sig_afs,length(Ch4_sig_positions)), get_geno, ne=2000)

ch1_out <- data.frame("chr"=rep("Ch1",Ch_length),"pos"=c(Ch1_neut_positions,Ch1_sig_positions),"ref"=rep('A',Ch_length),"alt"=rep('A/T',Ch_length),"geno"=c(Ch1_geno,Ch1_sig_geno))
ch2_out <- data.frame("chr"=rep("Ch2",Ch_length),"pos"=c(Ch2_neut_positions,Ch2_sig_positions),"ref"=rep('A',Ch_length),"alt"=rep('A/T',Ch_length),"geno"=c(Ch2_geno,Ch2_sig_geno))
ch3_out <- data.frame("chr"=rep("Ch3",Ch_length),"pos"=c(Ch3_neut_positions,Ch3_sig_positions),"ref"=rep('A',Ch_length),"alt"=rep('A/T',Ch_length),"geno"=c(Ch3_geno,Ch3_sig_geno))
ch4_out <- data.frame("chr"=rep("Ch4",Ch_length),"pos"=c(Ch4_neut_positions,Ch4_sig_positions),"ref"=rep('A',Ch_length),"alt"=rep('A/T',Ch_length),"geno"=c(Ch4_geno,Ch4_sig_geno))

hapout <- rbind(ch1_out,ch2_out,ch3_out,ch4_out)

write.table(hapout,"haplotypes.500Mb.ne2000.mimhap",sep='\t',quote=F,row.names=F,col.names = F)

# prep selected SNP selection coefficients
#sig <- filter(all, CMH_qval<0.05)
sig <- readRDS('lab.res.new.cmh05.lrt05.RDS')
selcoef <- sig$selCoef[which(sig$T0_AF != 0 & sig$T0_AF != 1)]

sel_out <- data.frame(c(rep("Ch1",length(Ch1_sig_positions)),rep("Ch2",length(Ch2_sig_positions)),rep("Ch3",length(Ch3_sig_positions)),rep("Ch4",length(Ch4_sig_positions))),
                        c(Ch1_sig_positions,Ch2_sig_positions,Ch3_sig_positions,Ch4_sig_positions),
                        rep("T/A",length(selcoef)-2),
                        selcoef[1:1136],
                        rep(0.5,length(selcoef)-2))

write.table(sel_out,"selection_coefficients.cmh05.lrt05.txt",sep='\t',quote=F,row.names=F,col.names = F)

##########
# one 10mb chromosome with 25 selected SNPs at real frequencies and selection selection_coefficients
library(dplyr)
all <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.new.RDS')
afs <- all$T0_AF_rising[which(all$T0_AF_rising != 0 & all$T0_AF_rising != 1)]

sig <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.new.cmh05.lrt05.RDS')
sig_afs <- sig$T0_AF_rising[which(sig$T0_AF_rising != 0 & sig$T0_AF_rising != 1)]
#set effective pop size
ne <- 250

get_geno <- function(af,ne){
  genos <- rbinom(2*ne,1,af)
  genos[which(genos==0)] <- 'T'
  genos[which(genos==1)] <- 'A'
  genos_h1 <- genos[1:ne]
  genos_h2 <- genos[(ne+1):(2*ne)]
  hap <- paste(genos_h1,genos_h2,sep='',collapse=' ')
  return(hap)
}

get_geno_singleton <- function(ne){
  hap1 <- rep("T",ne)
  hap2 <- rep("T",ne)
  hap1[sample(1:ne,1)] <- "A"
  hap <- paste(hap1,hap2,sep='',collapse=' ')
  return(hap)
}

Ch_length <- 7000
Ch1_positions <- sample(1:9999999, Ch_length)
Ch1_sig_positions <- sample(Ch1_positions,25)
Ch1_neut_positions <- Ch1_positions[which(!(Ch1_positions %in% Ch1_sig_positions))]
Ch1_geno <- sapply(sample(afs,length(Ch1_neut_positions)), get_geno, ne=ne)
Ch1_sig_geno <- sapply(sample(sig_afs,length(Ch1_sig_positions)), get_geno, ne=ne)
#Ch1_sig_geno <- sapply(rep(ne,25), get_geno_singleton)


hap_out <- data.frame("chr"=rep("Ch1",Ch_length),"pos"=c(Ch1_neut_positions,Ch1_sig_positions),"ref"=rep('A',Ch_length),"alt"=rep('A/T',Ch_length),"geno"=c(Ch1_geno,Ch1_sig_geno))
write.table(hap_out,"haplotypes.10Mb.ne250.rep10.mimhap",sep='\t',quote=F,row.names=F,col.names = F)

selcoef <- round(abs(sig$selCoef[which(sig$T0_AF_rising != 0 & sig$T0_AF_rising != 1)]),3)

sel_out <- data.frame(rep("Ch1",length(Ch1_sig_positions)),
                        Ch1_sig_positions,
                        rep("T/A",25),
                        sample(selcoef,25),
                        rep(0.5,25))

write.table(sel_out,"selection_coefficients.ne250.rep10.txt",sep='\t',quote=F,row.names=F,col.names = F)

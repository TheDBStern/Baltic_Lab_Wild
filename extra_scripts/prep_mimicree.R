##########
# one 10mb chromosome with 25 selected SNPs at real frequencies and selection selection_coefficients
##########
library(dplyr)
all <- readRDS('lab.all.RDS')
afs <- all$T0_AF_rising[which(all$T0_AF_rising != 0 & all$T0_AF_rising != 1)]

sig <- readRDS('lab.sig.RDS')
sig_afs <- sig$T0_AF_rising[which(sig$T0_AF_rising != 0 & sig$T0_AF_rising != 1)]

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

write.table(sel_out,"selection_coefficients.ne2000.rep1.txt",sep='\t',quote=F,row.names=F,col.names = F)

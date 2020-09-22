#Get starting frequencies of fixed SNPs

library(ggplot2)
library(dplyr)
dat <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/full_data.rawFreqs.RDS')
bse0_tmp <- filter(dat,Beaker=="BSE-0")
bse0_tmp$Treat <- rep("Treatment",nrow(bse0_tmp))
bse0b_tmp <- filter(dat,Beaker=="BSE-0B")
bse0b_tmp$Treat <- rep("Treatment",nrow(bse0b_tmp))
dat <- rbind(dat,bse0_tmp)
dat <- rbind(dat,bse0b_tmp)
dat$folded <- sapply(dat$value,fold_AF)


res <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.new.cmh05.RDS')
sig <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.new.cmh05.lrt05.RDS')
colnames(sig)[1:2] <- c('Transcript','Position')
sig_dat <- merge(dat,sig,by=c('Transcript','Position'))

## polarize so show SNP rising in freq in selected lines
sig_dat <- as.data.frame(sig_dat)

snpnums <- unique(sig_dat$SNP)

polarized <- data.frame()
for (i in 1:length(snpnums)){
	snpnum <- snpnums[i]
	print(snpnum)
	dt <- filter(sig_dat,SNP==snpnum)
	#mod <- lm(value~Generation,data=filter(dt,Treat=="Treatment"))
	slope <- dt$estS_lm[[1]]
	if (slope < 0 ){
		dt$value <- (1- dt$value)
		}
	polarized <- rbind(polarized, dt)
	}

dat <- filter(polarized,Treat=="Treatment" & Generation==10)

fixed <- filter(dat,value==1)
# get linear slope estimates for sig SNPs in selected lines

library(dplyr)
library(lme4)
library(parallel)

## for transformed divergence from ancestor
dat <- readRDS('full_data_lmm.RDS')
snpdet <- read.table('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.snpdet')
dat$Transcript <- rep(snpdet[,1],26)
dat$Position <- rep(snpdet[,2],26)

res <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.cmh05.RDS')
sig <- res[which(res$LRT_qval <= 0.05),]

sig_dat <- data.frame()

for (i in 1:nrow(sig)){
	trans <- toString(sig[i,1])
	pos <- sig[i,2]
	snp <- dat[which(dat$Transcript == trans & dat$Position == pos),]
	sig_dat <- rbind(sig_dat,snp)
	}

sig_snp_nums <- unique(sig_dat$SNP)

lm_test <-function(i){
  snpnum <- sig_snp_nums[i]
  print(snpnum)
  SNPdat= sig_dat[which(sig_dat$SNP == snpnum & Treat == "Treatment"),]
  neff <- (100*SNPdat$Coverage-1)/ (100+SNPdat$Coverage)
  mod=lmer(value~Generation+(1|Beaker), weights=neff, data=SNPdat,REML=F)
  sum <- summary(mod)
  return(coef(sum)[2,1])
}

lm_res <- lapply(1:length(sig_snp_nums),lm_test)

lm_res <- unlist(lm_res)

saveRDS(lm_res, 'sig_snps.transformed_divergence_correlation_coef.RDS')

## for raw freq frequencies

dat <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/full_data.rawFreqs.RDS')
bse0_tmp <- filter(dat,Beaker=="BSE-0")
bse0_tmp$Treat <- rep("Treatment",nrow(bse0_tmp))
bse0b_tmp <- filter(dat,Beaker=="BSE-0B")
bse0b_tmp$Treat <- rep("Treatment",nrow(bse0b_tmp))
dat <- rbind(dat,bse0_tmp)
dat <- rbind(dat,bse0b_tmp)
dat[which(dat$value==0),3] <- 0.01
dat[which(dat$value==1),3] <- 0.99

dat$AF_logit <- log(dat$value/(1-dat$value))
dat$AF_logit10 <- log10(dat$value/(1-dat$value))

#dat$AF_logit <- logit(dat$value,adjust=0.01)

res <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.cmh05.RDS')
sig <- res[which(res$LRT_qval <= 0.05),]

sig_dat <- data.frame()

for (i in 1:nrow(sig)){
	trans <- toString(sig[i,1])
	pos <- sig[i,2]
	snp <- dat[which(dat$Transcript == trans & dat$Position == pos),]
	sig_dat <- rbind(sig_dat,snp)
	}

sig_snp_nums <- unique(dat$SNP)

lm_test <-function(i){
  print(i)
  SNPdat= dat[which(dat$SNP == i & dat$Treat == "Treatment"),]
  neff <- (100*SNPdat$Coverage-1)/ (100+SNPdat$Coverage)
  mod=lmer(value~Generation+(1|Beaker), weights=neff, data=SNPdat,REML=F)
  sum <- summary(mod)
  return(coef(sum)[2,1])
}

value_lm_res <- lapply(1:length(sig_snp_nums),lm_test)

#logit_lm_res <- lapply(1:length(sig_snp_nums),lm_test)


lm_res <- unlist(lm_res)

saveRDS(lm_res, 'sig_snps.rawAF_correlation_coef.RDS')



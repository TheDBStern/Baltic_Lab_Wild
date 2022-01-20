library(data.table)
library(lme4)
library(parallel)
library(dplyr)


dat <- readRDS('hap_blocks.rawFreqs.RDS') ## generated with prep_lmm.rawFreqs.R script
## create duplicate of starting population for treatment lines
bse0_tmp <- filter(dat,Beaker=="BSE-0")
bse0_tmp$Treat <- rep("Treatment",nrow(bse0_tmp))
bse0b_tmp <- filter(dat,Beaker=="BSE-0B")
bse0b_tmp$Treat <- rep("Treatment",nrow(bse0b_tmp))
dat <- rbind(dat,bse0_tmp)
dat <- rbind(dat,bse0b_tmp)
## convert allele frequencies of 0 to 0.01
dat[which(dat$value==0),3] <- 0.01
dat[which(dat$value==1),3] <- 0.99

dat$AF_logit <- log(dat$value/(1-dat$value))

dat <- filter(dat,dat$Treat == "Treatment")

loc_nums <- unique(dat$loc)

estS <-function(i){
  print(i)
  locdat= dat[which(dat$loc == i),]
  neff <- (100*locdat$Coverage-1)/ (100+locdat$Coverage)
  mod=lmer(AF_logit~Generation+(1|Beaker), weights=neff, data=locdat,REML=F)
  sum <- summary(mod)
  return(coef(sum)[2,1])
}


res <- mclapply(1:length(loc_nums),estS,mc.cores=4)
res <- unlist(res)
saveRDS(res, 'estS.logit.lmm.RDS')

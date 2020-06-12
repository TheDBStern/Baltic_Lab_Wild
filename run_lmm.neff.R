library(data.table)
library(lme4)
library(parallel)
library(car)

dt <- readRDS('full_data_lmm.RDS')

lrt_test <-function(i){
  print(i)
  SNPdat= dt[which(dt$SNP == i),]
  neff <- (100*SNPdat$Coverage-1)/ (100+SNPdat$Coverage)
  mod1=lmer(value~Generation+(1|Beaker), weights=neff,data=SNPdat,REML=F)
  mod2=lmer(value~Treat*Generation+(1|Beaker), weights=neff, data=SNPdat,REML=F)
  aov = anova(mod1,mod2)
  return(aov$'Pr(>Chisq)'[[2]][1])
}

lrt_res <- mclapply(1:max(dt$SNP),lrt_test,mc.cores=8)

lrt_res <- unlist(lrt_res)

saveRDS(lrt_res, 'trtxgen.ranBeaker.neff.lrt.pvals.RDS')

## basic treat test

lm_test <-function(i){
  print(i)
  SNPdat= dt[which(dt$SNP == i),]
  neff <- (100*SNPdat$Coverage-1)/ (100+SNPdat$Coverage)
  mod=lm(value~Treat, weights=neff,data=SNPdat)
  sum <- summary(mod)
  return(sum$coefficients[2,4])
}

lm_res <- mclapply(1:max(dt$SNP),lm_test,mc.cores=8)

lm_res <- unlist(lm_res)

saveRDS(lm_res, 'basic_lm.trans.neff.pvals.RDS')


## levene test

lv_test <-function(i){
  print(i)
  SNPdat= dt[which(dt$SNP == i),]
  neff <- (100*SNPdat$Coverage-1)/ (100+SNPdat$Coverage)
  mod=leveneTest(value~Treat, weights=neff,data=SNPdat)
  return(mod$'Pr(>F)'[1])
}

lv_res <- mclapply(1:max(dt$SNP),lv_test,mc.cores=2)

lv_res <- unlist(lv_res)

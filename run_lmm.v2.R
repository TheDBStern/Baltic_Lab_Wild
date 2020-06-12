library(data.table)
library(lme4)
library(parallel)

dt <- readRDS('full_data_lmm.RDS')

lrt_test <-function(i){
  print(i)
  SNPdat= dt[which(dt$SNP == i),]
  mod1=lmer(value~Generation+Coverage+(1|Beaker), data=SNPdat,REML=F)
  mod2=lmer(value~Treat*Generation+Coverage+(1|Beaker), data=SNPdat,REML=F)
  aov = anova(mod1,mod2)
  return(aov$'Pr(>Chisq)'[[2]][1])
}

lrt_res <- mclapply(1:max(dt$SNP),lrt_test,mc.cores=8)

lrt_res <- unlist(lrt_res)

saveRDS(lrt_res, 'trtxgen+cov.ranBeaker.lrt.pvals.RDS')
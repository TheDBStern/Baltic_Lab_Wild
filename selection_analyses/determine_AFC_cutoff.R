library(poolSeq)


#neutral evolution for haplotype blocks
dat <- readRDS('data/hap_blocks.res.RDS')
res <- data.frame()
for (i in 1:nrow(dat)){
  dt <- dat[i,]
  sims <- wf.traj(p0=rep(dt$T0_AF_rising,10000), s=0, Ne=1750, t=c(0, 6,10))
  gen10_afc <- abs(sims[,3]- sims[,1])
  gen10_999 <- quantile(gen10_afc,0.999)

  df <- data.frame("tag"=dt$tag,"T0"=dt$T0_AF,
                          "cutoff"=gen10_999
  res <- rbind(res,df)
}

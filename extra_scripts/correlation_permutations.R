library(dplyr)

res <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.RDS')

cmh_p <- c()
cmh_cor <- c()
lrt_p <- c()
lrt_cor <- c()

for (i in 1:1000){
	print(i)
	dat <- sample_n(res,1156)
	cmh_res <- cor.test(dat$CMH_sel, dat$T0_AF,method="spearman")
	lrt_res <- cor.test(dat$LRT_pval, dat$T0_AF,method="spearman")
	cmh_p <- c(cmh_p, cmh_res$p.value)
	cmh_cor <- c(cmh_cor, cmh_res$estimate)
	lrt_p <- c(lrt_p, lrt_res$p.value)
	lrt_cor <- c(lrt_cor, lrt_res$estimate)
}

mean(cmh_p)
mean(cmh_cor)
mean(lrt_p)
mean(lrt_cor)

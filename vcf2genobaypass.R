library(poolfstat)

popnames = c("BSE-0","BSE-0B",
			 "BSE3C","BSE4C",
			"BSE-1-T1","BSE-1-T2","BSE-1-T2B",
			"BSE-2-T1","BSE-2-T2","BSE-2-T2B",
			"BSE-3-T1","BSE-3-T2",
			"BSE-4-T1","BSE-4-T2",
			"BSE-5-T1","BSE-5-T2",
			"BSE-6-T1","BSE-6-T2",
			"BSE-7-T1",
			"BSE-8-T1","BSE-8-T2",
			"BSE-9-T1","BSE-9-T2",
			"BSE-10-T1",
			"BSE-11-T1","BSE-11-T2",
			"BSE-12-T1","BSE-12-T2")

dat <- vcf2pooldata(vcf.file="lab.vcf",poolsizes=c(200,rep(100,20)),
			poolnames=popnames,min.cov.per.pool = 10,max.cov.per.pool=200,min.maf=0.01,nlines.per.readblock=1000000)

saveRDS(dat,"pooldat.RDS")

pooldata2genobaypass(dat,writing.dir=getwd())

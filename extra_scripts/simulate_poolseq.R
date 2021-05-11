# simulate pool-seq data under different starting frequencies and selection coefficients
# to match those estimated in the real data

library(poolSeq)
library(dplyr)

numSnps <- 1156
poolsize <- 100
freqdat <- read.table('../../lab.genobaypass_freq',h=F)
covdat <- read.table('../../lab.genobaypass_cov',h=F)
snpdet <- read.table('../../lab.snpdet',h=F)[,1:2]
colnames(snpdet) <- c("Pseudo_Transcript","Pseudo_Position")
freqdat <- cbind(freqdat,snpdet)
res_all <- readRDS('../../lab.res.new.RDS')
seldat <- readRDS('../../lab.res.new.cmh05.lrt05.RDS')
#non_seldat <- filter(res_all, LRT_pval > max(seldat$LRT_pval) | CMH_qval > 0.05)
selfreq <- merge(freqdat,seldat,by=c("Pseudo_Transcript","Pseudo_Position"))
allfreq <- merge(freqdat,res_all,by=c("Pseudo_Transcript","Pseudo_Position"))

Ne <- 2000


#########
## simulate 1156 loci under selection in 10 populations for ten replicate simulations
#########

for (rep in 11:100){
print(paste("Rep ",rep,sep=''))
res_freq <- c()
res_cov <- c()

for (i in 1:1156){
	#print(i)
	start <- rowMeans(sample_n(allfreq,1)[,3:4])
	while (start == 0) {start <- rowMeans(sample_n(allfreq,1)[,3:4]) }
	sel <- sample_n(seldat,1)$selCoef
	snp_af <- c()
	snp_cov <- c()
	#selected
	for (x in 1:10){
		simTraj <- wf.traj(p0=start, s=sel, Ne=Ne, t=c(0, 6, 10))
		simTraj <- sample.alleles(simTraj, size=poolsize, mode="individuals", Ncensus=500)
		af <- sample.alleles(simTraj, size=round(mean(as.numeric(sample_n(covdat,1)[,1:2]))), mode="coverage")
		simTraj <- af$p.smpld
		simCov <- af$size
		snp_af <- c(snp_af,simTraj)
		snp_cov <- c(snp_cov,simCov)
		}
	res_freq <- rbind(res_freq,snp_af)
	res_cov <- rbind(res_cov,snp_cov)
	}

write.table(res_freq,paste('emp_vals.1156sel.allStart.rep',rep,'.freq',sep=''),sep='\t',row.names=F,quote=F)
write.table(res_cov,paste('emp_vals.1156sel.allStart.rep',rep,'.cov',sep=''),sep='\t',row.names=F,quote=F)
}

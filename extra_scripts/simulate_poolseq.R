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
non_seldat <- filter(res_all, LRT_pval > max(seldat$LRT_pval) | CMH_qval > 0.05)
selfreq <- merge(freqdat,seldat,by=c("Pseudo_Transcript","Pseudo_Position"))

#nrep <- 10
Ne <- 2000

######### 
###all the same
########
res_freq <- c()
res_cov <- c()
for (i in 1:numSnps){
	print(i)
	#start <- runif(1, min = 0.01, max = 0.99)
	start <- rowMeans(sample_n(freqdat,1)[,1:2])
	#sel <- rnorm(1,mean = 0.119, sd = 0.0456)
	#sel <- 0
	sel <- sample_n(seldat,1)$estS_lm
	snp_af <- c()
	snp_cov <- c()
	for (x in 1:nrep){	
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
	
write.table(res_freq,'sim.freq',sep='\t',row.names=F,quote=F)
write.table(res_cov,'sim.cov',sep='\t',row.names=F,quote=F)

####### 
####control and sel
######
res_freq <- c()
res_cov <- c()

#simulate 10000 selected snps and 90000 non-selected snps
for (i in 1:1000){
	print(i)
	#start <- runif(1, min = 0.01, max = 0.99)
	start <- rowMeans(sample_n(freqdat,1)[,1:2])
	sel <- 0.119
	#sel <- rnorm(1,mean = 0.119, sd = 0.0456)
	#sel <- sample_n(seldat,1)$estS_lm
	#nonsel <- sample_n(non_seldat,1)$estS_lm
	nonsel <- 0
	snp_af <- c()
	snp_cov <- c()
	#control (i.e. neutral)
	for (x in 1:4){	
		simTraj <- wf.traj(p0=start, s=nonsel, Ne=Ne, t=c(0, 6, 10))
		simTraj <- sample.alleles(simTraj, size=poolsize, mode="individuals", Ncensus=500)
		af <- sample.alleles(simTraj, size=round(mean(as.numeric(sample_n(covdat,1)[,1:2]))), mode="coverage")
		simTraj <- af$p.smpld
		simCov <- af$size
		snp_af <- c(snp_af,simTraj)
		snp_cov <- c(snp_cov,simCov)
		}
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
for (i in 1:9000){
	print(i)
	#start <- runif(1, min = 0.01, max = 0.99)      
	start <- rowMeans(sample_n(freqdat,1)[,1:2]) 
	nonsel <- 0
	#nonsel <- sample_n(non_seldat,1)$estS_lm
	snp_af <- c()
	snp_cov <- c()
	#all neutral
	for (x in 1:14){	
		simTraj <- wf.traj(p0=start, s=nonsel, Ne=Ne, t=c(0, 6, 10))
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


write.table(res_freq,'sim_4ctrl_10sel.1ksel_9knon.freq',sep='\t',row.names=F,quote=F)
write.table(res_cov,'sim_4ctrl_10sel.1ksel_9knon.cov',sep='\t',row.names=F,quote=F)

########
## test how small frequency changes can be detected by our design
########

res_freq <- c()
res_cov <- c()

#simulate 100000 selected snps 
for (i in 1:1156){
	print(i)
	#start <- runif(1, min = 0.01, max = 0.99)
	start <- rowMeans(sample_n(selfreq,1)[,3:4])
	#start <- rowMeans(sample_n(freqdat,1)[,1:2])
	while (start == 0) {start <- rowMeans(sample_n(selfreq,1)[,3:4]) }
	#sel <- 0.119
	#sel <- rnorm(1,mean = 0.119, sd = 0.0456)
	sel <- sample_n(seldat,1)$estS_lm
	#sel <- runif(1,min=0.01,max=0.5)
	#nonsel <- sample_n(non_seldat,1)$estS_lm
	nonsel <- 0
	snp_af <- c()
	snp_cov <- c()
	#control (i.e. neutral)
	for (x in 1:4){	
		simTraj <- wf.traj(p0=start, s=nonsel, Ne=Ne, t=c(0, 6, 10))
		simTraj <- sample.alleles(simTraj, size=poolsize, mode="individuals", Ncensus=500)
		af <- sample.alleles(simTraj, size=round(mean(as.numeric(sample_n(covdat,1)[,1:2]))), mode="coverage")
		simTraj <- af$p.smpld
		simCov <- af$size
		snp_af <- c(snp_af,simTraj)
		snp_cov <- c(snp_cov,simCov)
		}
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
for (i in 1:1000){
	print(i)
	#start <- runif(1, min = 0.01, max = 0.99)      
	start <- rowMeans(sample_n(selfreq,1)[,3:4])
	while (start == 0) {start <- rowMeans(sample_n(selfreq,1)[,3:4]) }
	nonsel <- 0
	#nonsel <- sample_n(non_seldat,1)$estS_lm
	snp_af <- c()
	snp_cov <- c()
	#all neutral
	for (x in 1:14){	
		simTraj <- wf.traj(p0=start, s=nonsel, Ne=Ne, t=c(0, 6, 10))
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

write.table(res_freq,'emp_vals.1156sel_1knon.selStart.freq',sep='\t',row.names=F,quote=F)
write.table(res_cov,'emp_vals.1156sel_1knon.selStart.cov',sep='\t',row.names=F,quote=F)


#########
## simulate 1156 loci under selection in 10 populations for ten replicate simulations
#########

for (rep in 1:10){
print(paste("Rep ",rep,sep=''))
res_freq <- c()
res_cov <- c()

for (i in 1:1156){
	#print(i)
	start <- rowMeans(sample_n(selfreq,1)[,3:4])
	while (start == 0) {start <- rowMeans(sample_n(selfreq,1)[,3:4]) }
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

write.table(res_freq,paste('emp_vals.1156sel.selStart.rep',rep,'.freq',sep=''),sep='\t',row.names=F,quote=F)
write.table(res_cov,paste('emp_vals.1156sel.selStart.rep',rep,'.cov',sep=''),sep='\t',row.names=F,quote=F)
}

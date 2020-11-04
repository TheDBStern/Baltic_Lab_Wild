library(ACER)
library(matrixStats)
freqs <- read.table("../lab.genobaypass_freq")
covs <- read.table("../lab.genobaypass_cov")

########### CMH test #####################
###############
### T0,T1, T2
###############
freqs$T0 <- rowMeans(freqs[,1:2])
covs$T0 <- rowSums(covs[,1:2])

#selected
freqdat <- as.matrix(freqs[,c(29,11,12,29,13,14,29,15,16,29,17,18,29,20,21,29,22,23,29,25,26,29,27,28)])
covdat <- as.matrix(covs[,c(29,11,12,29,13,14,29,15,16,29,17,18,29,20,21,29,22,23,29,25,26,29,27,28)])
rep <- c(1:10)
ps <- rep(c(200,100,100), length(rep))
tp <- c(0,6,10)

res <- adapted.cmh.test(freq=freqdat, coverage=covdat,IntGen=TRUE, order=0, gen=tp,Ne=rep(2000,10), repl=rep, poolSize=ps,RetVal = 2)
saveRDS(res, "acer.cmh.Ne100.sel.RDS")

#controls
freqs$BSE1_T2 <- rowMeans(freqs[,6:7])
covs$BSE1_T2 <- rowSums(covs[,6:7])
freqs$BSE2_T2 <- rowMeans(freqs[,9:10])
covs$BSE2_T2 <- rowSums(covs[,9:10])

freqdat <- as.matrix(freqs[,c(29,5,30,29,8,30)])
covdat <- as.matrix(covs[,c(29,5,30,29,8,30)])
rep <- c(1:2)
ps <- rep(c(200,100,200), length(rep))
tp <- c(0,6,10)

res_ctrl <- adapted.cmh.test(freq=freqdat, coverage=covdat,IntGen=TRUE, order=0, gen=tp,Ne=rep(1000,2), repl=rep, poolSize=ps,RetVal = 2)
saveRDS(res_ctrl, "acer.cmh.Ne1000.ctrls.RDS")

########### CMH test for simulated data #####################
###############
### T0,T1, T2
###############
freqs <- read.table("unifsel.10ksel_90knon.freq",h=T)
covs <- read.table("unifsel.10ksel_90knon.cov",h=T)
freqs$T0 <- rowMeans(freqs[,c(1,4)])
covs$T0 <- rowSums(covs[,c(1,4)])

#selected
freqdat <- as.matrix(freqs[,c(43,14,15,43,17,18,43,20,21,43,23,24,43,26,27,43,29,30,43,32,33,43,35,36)])
covdat <- as.matrix(covs[,c(43,14,15,43,17,18,43,20,21,43,23,24,43,26,27,43,29,30,43,32,33,43,35,36)])
rep <- c(1:8)
ps <- rep(c(200,100,100), length(rep))
tp <- c(0,6,10)

res <- adapted.cmh.test(freq=freqdat, coverage=covdat,IntGen=TRUE, order=0, gen=tp,Ne=rep(2000,8), repl=rep, poolSize=ps,RetVal = 2)

saveRDS(res, "cmh_res.RDS")

###############
### T1 vs T2, selected
###############
freqdat <- as.matrix(freqs[,c(11,12,13,14,15,16,17,18,20,21,22,23,25,26,27,28)])
covdat <- as.matrix(covs[,c(11,12,13,14,15,16,17,18,20,21,22,23,25,26,27,28)])
rep <- c(1:8)
ps <- rep(100, 2*length(rep))
tp <- c(6,10)

res <- adapted.cmh.test(freq=freqdat, coverage=covdat,IntGen=FALSE, order=0, gen=tp,Ne=rep(1000,8), repl=rep, poolSize=ps,RetVal = 2)
saveRDS(res, "acer.cmh.pool_correction_only.ctrls.RDS")



######## Chi-square test ###################
freqdat <- read.table('../lab.genobaypass_freq',h=F)
covdat <- read.table('../lab.genobaypass_cov',h=F)
tp <- c(0,6,10)
ps <- c(200,100,100)

BSE1_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,5],F10=rowMeans(cbind(freqdat[,6:7])))
BSE1_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,5],F10=round(rowSums(cbind(covdat[,6:7])),0))
BSE2_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,8],F10=rowMeans(cbind(freqdat[,9:10])))
BSE2_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,8],F10=round(rowSums(cbind(covdat[,9:10])),0))
BSE3_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,11],F10=freqdat[,12])
BSE3_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,11],F10=covdat[,12])
BSE4_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,13],F10=freqdat[,14])
BSE4_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,13],F10=covdat[,14])
BSE5_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,15],F10=freqdat[,16])
BSE5_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,15],F10=covdat[,16])
BSE6_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,17],F10=freqdat[,18])
BSE6_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,17],F10=covdat[,18])
BSE8_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,20],F10=freqdat[,21])
BSE8_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,20],F10=covdat[,21])
BSE9_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,22],F10=freqdat[,23])
BSE9_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,22],F10=covdat[,23])
BSE11_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,25],F10=freqdat[,26])
BSE11_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,25],F10=covdat[,26])
BSE12_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,27],F10=freqdat[,28])
BSE12_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,27],F10=covdat[,28])

#Estimate Ne with LLS
BSE1_Ne <- estimateNe(p0=BSE1_freq[,"F0"], pt=BSE1_freq[,"F6"], cov0=BSE1_cov[,"F0"], covt=BSE1_cov[,"F6"], t=6, poolSize=c(200, 100), method="P.planII",truncAF=0.01)
BSE2_Ne <- estimateNe(p0=BSE2_freq[,"F0"], pt=BSE2_freq[,"F6"], cov0=BSE2_cov[,"F0"], covt=BSE2_cov[,"F6"], t=6, poolSize=c(200, 100), method="P.planII",truncAF=0.01)

BSE3_Ne <- estimateNe(p0=BSE3_freq[,"F0"], pt=BSE3_freq[,"F10"], cov0=BSE3_cov[,"F0"], covt=BSE3_cov[,"F10"], t=10, poolSize=c(200, 100), method="P.planI",truncAF=0.01, Ncensus=c(500,500))
BSE4_Ne <- estimateNe(p0=BSE4_freq[,"F0"], pt=BSE4_freq[,"F10"], cov0=BSE4_cov[,"F0"], covt=BSE4_cov[,"F10"], t=10, poolSize=c(200, 100), method="P.planI",truncAF=0.01,Ncensus=c(500,500))
BSE5_Ne <- estimateNe(p0=BSE5_freq[,"F0"], pt=BSE5_freq[,"F10"], cov0=BSE5_cov[,"F0"], covt=BSE5_cov[,"F10"], t=10, poolSize=c(200, 100), method="P.planI",truncAF=0.01,Ncensus=c(500,500))
BSE6_Ne <- estimateNe(p0=BSE6_freq[,"F0"], pt=BSE6_freq[,"F10"], cov0=BSE6_cov[,"F0"], covt=BSE6_cov[,"F10"], t=10, poolSize=c(200, 100), method="P.planI",truncAF=0.01,Ncensus=c(500,500))
BSE8_Ne <- estimateNe(p0=BSE8_freq[,"F0"], pt=BSE8_freq[,"F10"], cov0=BSE8_cov[,"F0"], covt=BSE8_cov[,"F10"], t=10, poolSize=c(200, 100), method="P.planI",truncAF=0.01,Ncensus=c(500,500))
BSE9_Ne <- estimateNe(p0=BSE9_freq[,"F0"], pt=BSE9_freq[,"F10"], cov0=BSE9_cov[,"F0"], covt=BSE9_cov[,"F10"], t=10, poolSize=c(200, 100), method="P.planI",truncAF=0.01,Ncensus=c(500,500))
BSE11_Ne <- estimateNe(p0=BSE11_freq[,"F0"], pt=BSE11_freq[,"F10"], cov0=BSE11_cov[,"F0"], covt=BSE11_cov[,"F10"], t=10, poolSize=c(200, 100), method="P.planI",truncAF=0.01,Ncensus=c(500,500))
BSE12_Ne <- estimateNe(p0=BSE12_freq[,"F0"], pt=BSE12_freq[,"F10"], cov0=BSE12_cov[,"F0"], covt=BSE12_cov[,"F10"], t=10, poolSize=c(200, 100), method="P.planI",truncAF=0.01,Ncensus=c(500,500))

# estimate snps under selection with adapted chi-square

BSE3_res <- adapted.chisq.test(freq=BSE3_freq, coverage=BSE3_cov,IntGen=TRUE, gen=tp,Ne=BSE3_Ne, poolSize=ps,RetVal = 2)
BSE4_res <- adapted.chisq.test(freq=BSE4_freq, coverage=BSE4_cov,IntGen=TRUE, gen=tp,Ne=BSE4_Ne, poolSize=ps,RetVal = 2)
BSE5_res <- adapted.chisq.test(freq=BSE5_freq, coverage=BSE5_cov,IntGen=TRUE, gen=tp,Ne=BSE5_Ne, poolSize=ps,RetVal = 2)
BSE6_res <- adapted.chisq.test(freq=BSE6_freq, coverage=BSE6_cov,IntGen=TRUE, gen=tp,Ne=BSE6_Ne, poolSize=ps,RetVal = 2)
BSE8_res <- adapted.chisq.test(freq=BSE8_freq, coverage=BSE8_cov,IntGen=TRUE, gen=tp,Ne=BSE8_Ne, poolSize=ps,RetVal = 2)
BSE9_res <- adapted.chisq.test(freq=BSE9_freq, coverage=BSE9_cov,IntGen=TRUE, gen=tp,Ne=BSE9_Ne, poolSize=ps,RetVal = 2)
BSE11_res <- adapted.chisq.test(freq=BSE11_freq, coverage=BSE11_cov,IntGen=TRUE, gen=tp,Ne=BSE11_Ne, poolSize=ps,RetVal = 2)
BSE12_res <- adapted.chisq.test(freq=BSE12_freq, coverage=BSE12_cov,IntGen=TRUE, gen=tp,Ne=BSE12_Ne, poolSize=ps,RetVal = 2)

res_all <- data.frame(-log10(BSE3_res[,2]), -log10(BSE4_res[,2]),-log10(BSE5_res[,2]), -log10(BSE6_res[,2]),
		 -log10(BSE8_res[,2]), -log10(BSE9_res[,2]), -log10(BSE11_res[,2]), -log10(BSE12_res[,2]))

res_all <- data.frame(BSE3_res[,2], BSE4_res[,2], BSE5_res[,2], BSE6_res[,2],BSE8_res[,2], BSE9_res[,2], BSE11_res[,2], BSE12_res[,2])


######### CMH test using mimicree simulations
library(ACER)
library(haploReconstruct)
library(dplyr)
library(ggplot2)

reps <- 5
ne <- 250
dat <- sync_to_frequencies(file="ne250.rep10.5reps.sync",base.pops=rep(c(TRUE,rep(FALSE,2)),times=reps), header=F)

freqs <- as.matrix(dat[,7:(3*reps+6)])
cov <- matrix(
   rep(ne*2,nrow(freqs)*3*reps),
  nrow=nrow(freqs),
  ncol=3*reps,
 byrow = TRUE)

res <- adapted.cmh.test(freq=freqs, coverage=cov,poolSize = NULL,IntGen=TRUE, order=0, gen=c(0,6,10), Ne=rep(ne,reps), repl=1:reps, RetVal = 2)
#plot(-log10(res[,2]))
saveRDS(res, 'rep10.cmh.RDS')

true_sel <- read.table('selection_coefficients.ne250.rep10.txt',skip=1)
dat$pval <- -log10(res[,2])
dat$true <- ifelse(dat$pos %in% true_sel$V2, "Sel","Neutral")

top25 <- dat %>%
							arrange(desc(pval)) %>%
							head(25)
tp <- nrow(filter(top25, true=="Sel")) / 25
fp <- nrow(filter(top25, true=="Neutral")) / 25

write.table(cbind(tp,fp),'rep10.tpr_fpr.txt',sep='\t', quote=F,row.names=F)


p <- ggplot(dat, aes(x=pos, y=pval, color=true))
	p + geom_point(alpha=0.8) +
	scale_color_manual(values=c("black", "red")) +
	theme_bw()

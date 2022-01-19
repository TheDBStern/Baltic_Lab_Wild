library(ACER)
library(matrixStats)
library(qvalue)
library(dplyr)

freqs <- read.table("data/lab.genobaypass_freq")
covs <- read.table("data/lab.genobaypass_cov")
snpdet <- read.table('data/lab.snpdet')
########### CMH test #####################
# collect T0 data
freqs$T0 <- rowMeans(freqs[,1:2])
covs$T0 <- rowSums(covs[,1:2])

#selection lines only
#only 8 replicate lines that survived to gen 10
freqdat <- as.matrix(freqs[,c(29,11,12,29,13,14,29,15,16,29,17,18,29,20,21,29,22,23,29,25,26,29,27,28)])
covdat <- as.matrix(covs[,c(29,11,12,29,13,14,29,15,16,29,17,18,29,20,21,29,22,23,29,25,26,29,27,28)])
rep <- c(1:8)
ps <- rep(c(200,100,100), length(rep)) #poolsize for gen0 is 200 because 2 samples of 50 individuals
tp <- c(0,6,10)

Nes <- c(2271, 1873,1881,1659,1920,1582,1512,1502) #estimated effective pop sizes
cmh_res <- adapted.cmh.test(freq=freqdat, coverage=covdat,IntGen=TRUE, order=0, gen=tp,Ne=Nes, repl=rep, poolSize=ps,RetVal = 2)
res <- data.frame(snpdet[,1],snpdet[,2],-log10(cmh_res[,2]))
colnames(res) <- c("chr","pos","score")
saveRDS(res, "acer.cmh.RDS")
res$qval <- qvalue(10^(-res$score))$qvalues
sig <- filter(res, qval < 0.05)
saveRDS(sig,"acer.cmh.qval05.RDS")

######## Chi-square tests ###################
freqdat <- freqs
covdat <- covs

tp <- c(0,6,10)
ps <- c(200,100,100)
Nes <- c(2271, 1873,1881,1659,1920,1582,1512,1502)

#note: for gen 0, take mean frequencies and sum of coverage for 2 samples from starting population
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

# determine snps under selection with adapted chi-square
BSE3_res <- adapted.chisq.test(freq=BSE3_freq, coverage=BSE3_cov,IntGen=TRUE, gen=tp,Ne=Nes[[1]], poolSize=ps,RetVal = 2)
BSE4_res <- adapted.chisq.test(freq=BSE4_freq, coverage=BSE4_cov,IntGen=TRUE, gen=tp,Ne=Nes[[2]], poolSize=ps,RetVal = 2)
BSE5_res <- adapted.chisq.test(freq=BSE5_freq, coverage=BSE5_cov,IntGen=TRUE, gen=tp,Ne=Nes[[3]], poolSize=ps,RetVal = 2)
BSE6_res <- adapted.chisq.test(freq=BSE6_freq, coverage=BSE6_cov,IntGen=TRUE, gen=tp,Ne=Nes[[4]], poolSize=ps,RetVal = 2)
BSE8_res <- adapted.chisq.test(freq=BSE8_freq, coverage=BSE8_cov,IntGen=TRUE, gen=tp,Ne=Nes[[5]], poolSize=ps,RetVal = 2)
BSE9_res <- adapted.chisq.test(freq=BSE9_freq, coverage=BSE9_cov,IntGen=TRUE, gen=tp,Ne=Nes[[6]], poolSize=ps,RetVal = 2)
BSE11_res <- adapted.chisq.test(freq=BSE11_freq, coverage=BSE11_cov,IntGen=TRUE, gen=tp,Ne=Nes[[7]], poolSize=ps,RetVal = 2)
BSE12_res <- adapted.chisq.test(freq=BSE12_freq, coverage=BSE12_cov,IntGen=TRUE, gen=tp,Ne=Nes[[8]], poolSize=ps,RetVal = 2)

#calculate q values and compile results
res_all <- data.frame(snpdet[,1],snpdet[,2],
			qvalue(BSE3_res[,2])$qvalues, qvalue(BSE4_res[,2])$qvalues, qvalue(BSE5_res[,2])$qvalues,
			qvalue(BSE6_res[,2])$qvalues,qvalue(BSE8_res[,2])$qvalues, qvalue(BSE9_res[,2])$qvalues,
			qvalue(BSE11_res[,2])$qvalues, qvalue(BSE12_res[,2])$qvalues)
colnames(res_all) <- c("chr","pos","BSE3","BSE4","BSE5","BSE6","BSE8","BSE9","BSE11","BSE12")
sig <- res_all %>%
  filter_all(any_vars(. < 0.05))

saveRDS(sig, "acer.chisq.any_qval05.RDS")

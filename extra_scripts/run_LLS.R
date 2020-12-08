library(poolSeq)
library(parallel)

freqdat <- read.table('../lab.genobaypass_freq',h=F)
covdat <- read.table('../lab.genobaypass_cov',h=F)

### prep datasets for each replicate
BSE3_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,11],F10=freqdat[,12])
BSE4_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,13],F10=freqdat[,14])
BSE5_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,15],F10=freqdat[,16])
BSE6_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,17],F10=freqdat[,18])
BSE7_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,19])
BSE8_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,20],F10=freqdat[,21])
BSE9_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,22],F10=freqdat[,23])
BSE10_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,24])
BSE11_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,25],F10=freqdat[,26])
BSE12_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,27],F10=freqdat[,28])

BSE3_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,11],F10=covdat[,12])
BSE4_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,13],F10=covdat[,14])
BSE5_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,15],F10=covdat[,16])
BSE6_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,17],F10=covdat[,18])
BSE7_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,19])
BSE8_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,20],F10=covdat[,21])
BSE9_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,22],F10=covdat[,23])
BSE10_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F6=freqdat[,24])
BSE11_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,25],F10=covdat[,26])
BSE12_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F6=covdat[,27],F10=covdat[,28])

#s

s_res <- c()
for (i in 1:nrow(freqdat)){
	print(i)
	freq <- as.numeric(c(BSE3_freq[i,],BSE4_freq[i,],BSE5_freq[i,],BSE6_freq[i,],BSE8_freq[i,],BSE9_freq[i,],BSE11_freq[i,],
			BSE12_freq[i,],BSE7_freq[i,],BSE10_freq[i,]))
	res <- estimateSH(freq,Ne=2000,t=c(rep(c(0,6,10),8),rep(c(0,6),2)),h=0.5)
	s_res <- c(s_res,res$s)
	}

## attempt 2
freqdat <- read.table('../lab.genobaypass_freq',h=F)
covdat <- read.table('../lab.genobaypass_cov',h=F)
poporder <- c(1,2,11,13,15,17,19,20,22,24,25,27,12,14,16,18,21,23,26,28) #start (x2), gen6 (x10), gen10 (x8)

s_res <- c()
for (i in 1:nrow(freqdat)){
	print(i)
	freq <- as.numeric(freqdat[i,poporder])
	freq[which(freq==0)] <- 0.01
	freq[which(freq==1)] <- 0.99
	cov <- as.numeric(covdat[i,poporder])
	res <- estimateSH(freq,Ne=2000,t=c(0,0,rep(6,10),rep(10,8)),cov=cov,h=0.5, method="automatic")
	s_res <- c(s_res,res$s)
	}


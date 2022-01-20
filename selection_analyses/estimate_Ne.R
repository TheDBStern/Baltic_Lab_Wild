library(poolSeq)
library(parallel)
library(dplyr)

freqdat <- read.table('data/lab.snps_freq',h=F)
covdat <- read.table('data/lab.snps_cov',h=F)
snpdet <- read.table('data/lab.snpdet',h=F)

### prep datasets for each replicate
BS3C_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F20=freqdat[,3])
BS3C_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F20=covdat[,3])
BS4C_freq <- data.frame(F0=rowMeans(cbind(freqdat[,1:2])), F20=freqdat[,4])
BS4C_cov <- data.frame(F0=rowSums(cbind(covdat[,1:2])), F20=covdat[,4])
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

pops <- c("BS3C","BS4C","BSE1","BSE2","BSE3","BSE4","BSE5","BSE6","BSE8","BSE9","BSE11","BSE12")
cen <- c(500,1000,1500,2000,2500,3000,3500,4000,4500,5000)

### make estimate using all SNPs together
estimateNe(p0=BSE3_freq[,"F0"], pt=BSE3_freq[,"F10"], cov0=BSE3_cov[,"F0"], covt=BSE3_cov[,"F10"], t=10, poolSize=c(200, 100), method="P.planII",truncAF=0.05)

### make estimate in 10kb windows and taking median
for (i in cen){
win_ne <- estimateWndNe(snpdet[,1], snpdet[,2], 1000,
              p0=BSE3_freq[,"F0"], pt=BSE3_freq[,"F10"],
              cov0=BSE3_cov[,"F0"], covt=BSE3_cov[,"F10"], t=10,
              unit = "SNP", truncAF = 0.05,
              method = "P.planI", Ncensus=i,
              poolSize=c(200, 100))
print(median(filter(win_ne,Np.planI>0)$Np.planI))
}

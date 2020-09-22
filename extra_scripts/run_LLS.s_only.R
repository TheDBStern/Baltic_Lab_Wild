library(poolSeq)
library(parallel)

freqdat <- read.table('../lab.genobaypass_freq',h=F)
covdat <- read.table('../lab.genobaypass_cov',h=F)

### prep datasets for each replicate
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


## estimate Ne and s
#Ne
BSE1_Ne <- estimateNe(p0=BSE1_freq[,"F0"], pt=BSE1_freq[,"F6"], cov0=BSE1_cov[,"F0"], covt=BSE1_cov[,"F6"], t=6, poolSize=c(200, 100), Ncensus = 500, method="P.planI",truncAF=0.01)
BSE2_Ne <- estimateNe(p0=BSE2_freq[,"F0"], pt=BSE2_freq[,"F6"], cov0=BSE2_cov[,"F0"], covt=BSE2_cov[,"F6"], t=6, poolSize=c(200, 100), Ncensus = 500, method="P.planI",truncAF=0.01)
BSE3_Ne <- estimateNe(p0=BSE3_freq[,"F0"], pt=BSE3_freq[,"F6"], cov0=BSE3_cov[,"F0"], covt=BSE3_cov[,"F6"], t=6, poolSize=c(200, 100), Ncensus = 500, method="P.planI",truncAF=0.01)
BSE4_Ne <- estimateNe(p0=BSE4_freq[,"F0"], pt=BSE4_freq[,"F6"], cov0=BSE4_cov[,"F0"], covt=BSE4_cov[,"F6"], t=6, poolSize=c(200, 100), Ncensus = 500, method="P.planI",truncAF=0.01)
BSE5_Ne <- estimateNe(p0=BSE5_freq[,"F0"], pt=BSE5_freq[,"F6"], cov0=BSE5_cov[,"F0"], covt=BSE5_cov[,"F6"], t=6, poolSize=c(200, 100), Ncensus = 500, method="P.planI",truncAF=0.01)
BSE6_Ne <- estimateNe(p0=BSE6_freq[,"F0"], pt=BSE6_freq[,"F6"], cov0=BSE6_cov[,"F0"], covt=BSE6_cov[,"F6"], t=6, poolSize=c(200, 100), Ncensus = 500, method="P.planI",truncAF=0.01)
BSE8_Ne <- estimateNe(p0=BSE8_freq[,"F0"], pt=BSE8_freq[,"F6"], cov0=BSE8_cov[,"F0"], covt=BSE8_cov[,"F6"], t=6, poolSize=c(200, 100), Ncensus = 500, method="P.planI",truncAF=0.01)
BSE9_Ne <- estimateNe(p0=BSE9_freq[,"F0"], pt=BSE9_freq[,"F6"], cov0=BSE9_cov[,"F0"], covt=BSE9_cov[,"F6"], t=6, poolSize=c(200, 100), Ncensus = 500, method="P.planI",truncAF=0.01)
BSE11_Ne <- estimateNe(p0=BSE11_freq[,"F0"], pt=BSE11_freq[,"F6"], cov0=BSE11_cov[,"F0"], covt=BSE11_cov[,"F6"], t=6, poolSize=c(200, 100), Ncensus = 500, method="P.planI",truncAF=0.01)
BSE12_Ne <- estimateNe(p0=BSE12_freq[,"F0"], pt=BSE12_freq[,"F6"], cov0=BSE12_cov[,"F0"], covt=BSE12_cov[,"F6"], t=6, poolSize=c(200, 100), Ncensus = 500, method="P.planI",truncAF=0.01)

#s
run_LLS <- function(i,freqdat,covdat,Ne){
	freq <- freqdat[i,]
	cov <- covdat[i,]
	#res <- estimateSH(freq,Ne=abs(round(Ne,0)),t=c(0,6,10),h=0.5,cov=cov, N.ctraj = 1000,haploid=FALSE,simulate.p.value=T)
	res <- estimateSH(freq,Ne=abs(round(Ne,0)),t=c(0,6,10),h=0.5,cov=cov)

	print(i)
	#return(c(res$s,res$p.value))
	return(res$s)
	}


BSE1_res <- mclapply(1:nrow(freqdat),run_LLS,BSE1_freq,BSE1_cov,round(abs(BSE1_Ne),0),mc.cores=4)
BSE1_res <- data.frame(matrix(unlist(BSE1_res), nrow=nrow(freqdat), byrow=T),stringsAsFactors=FALSE)
colnames(BSE1_res) <- c("s")
write.table(BSE1_res,"BSE1.PplanI.Ncen500.txt",sep='\t',row.names=F,quote=F)

BSE2_res <- mclapply(1:nrow(freqdat),run_LLS,BSE2_freq,BSE2_cov,round(abs(BSE2_Ne),0),mc.cores=4)
BSE2_res <- data.frame(matrix(unlist(BSE2_res), nrow=nrow(freqdat), byrow=T),stringsAsFactors=FALSE)
colnames(BSE2_res) <- c("s")
write.table(BSE2_res,"BSE2.PplanI.Ncen500.txt",sep='\t',row.names=F,quote=F)

BSE3_res <- mclapply(1:nrow(freqdat),run_LLS,BSE3_freq,BSE3_cov,round(abs(BSE3_Ne),0),mc.cores=4)
BSE3_res <- data.frame(matrix(unlist(BSE3_res), nrow=nrow(freqdat), byrow=T),stringsAsFactors=FALSE)
colnames(BSE3_res) <- c("s")
write.table(BSE3_res,"BSE3.PplanI.Ncen500.txt",sep='\t',row.names=F,quote=F)

BSE4_res <- mclapply(1:nrow(freqdat),run_LLS,BSE4_freq,BSE4_cov,round(abs(BSE4_Ne),0),mc.cores=4)
BSE4_res <- data.frame(matrix(unlist(BSE4_res), nrow=nrow(freqdat), byrow=T),stringsAsFactors=FALSE)
colnames(BSE4_res) <- c("s")
write.table(BSE4_res,"BSE4.PplanI.Ncen500.txt",sep='\t',row.names=F,quote=F)

BSE5_res <- mclapply(1:nrow(freqdat),run_LLS,BSE5_freq,BSE5_cov,round(abs(BSE5_Ne),0),mc.cores=4)
BSE5_res <- data.frame(matrix(unlist(BSE5_res), nrow=nrow(freqdat), byrow=T),stringsAsFactors=FALSE)
colnames(BSE5_res) <- c("s")
write.table(BSE5_res,"BSE5.PplanI.Ncen500.txt",sep='\t',row.names=F,quote=F)

BSE6_res <- mclapply(1:nrow(freqdat),run_LLS,BSE6_freq,BSE6_cov,round(abs(BSE6_Ne),0),mc.cores=4)
BSE6_res <- data.frame(matrix(unlist(BSE6_res), nrow=nrow(freqdat), byrow=T),stringsAsFactors=FALSE)
colnames(BSE6_res) <- c("s")
write.table(BSE6_res,"BSE6.PplanI.Ncen500.txt",sep='\t',row.names=F,quote=F)

BSE8_res <- mclapply(1:nrow(freqdat),run_LLS,BSE8_freq,BSE8_cov,round(abs(BSE8_Ne),0),mc.cores=4)
BSE8_res <- data.frame(matrix(unlist(BSE8_res), nrow=nrow(freqdat), byrow=T),stringsAsFactors=FALSE)
colnames(BSE8_res) <- c("s")
write.table(BSE8_res,"BSE8.PplanI.Ncen500.txt",sep='\t',row.names=F,quote=F)

BSE9_res <- mclapply(1:nrow(freqdat),run_LLS,BSE9_freq,BSE9_cov,round(abs(BSE9_Ne),0),mc.cores=4)
BSE9_res <- data.frame(matrix(unlist(BSE9_res), nrow=nrow(freqdat), byrow=T),stringsAsFactors=FALSE)
colnames(BSE9_res) <- c("s")
write.table(BSE9_res,"BSE9.PplanI.Ncen500.txt",sep='\t',row.names=F,quote=F)

BSE11_res <- mclapply(1:nrow(freqdat),run_LLS,BSE11_freq,BSE11_cov,round(abs(BSE11_Ne),0),mc.cores=4)
BSE11_res <- data.frame(matrix(unlist(BSE11_res), nrow=nrow(freqdat), byrow=T),stringsAsFactors=FALSE)
colnames(BSE11_res) <- c("s")
write.table(BSE11_res,"BSE11.PplanI.Ncen500.txt",sep='\t',row.names=F,quote=F)

BSE12_res <- mclapply(1:nrow(freqdat),run_LLS,BSE12_freq,BSE12_cov,round(abs(BSE12_Ne),0),mc.cores=4)
BSE12_res <- data.frame(matrix(unlist(BSE12_res), nrow=nrow(freqdat), byrow=T),stringsAsFactors=FALSE)
colnames(BSE12_res) <- c("s")
write.table(BSE12_res,"BSE12.PplanI.Ncen500.txt",sep='\t',row.names=F,quote=F)



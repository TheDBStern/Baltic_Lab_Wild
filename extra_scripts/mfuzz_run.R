library(Mfuzz)
library(dplyr)

# prep data file
sig <- readRDS('../lab.res.new.cmh05.lrt05.RDS')[,c(1,2,24)]
sig <- cbind("GENE.ID"=paste(sig[,1],sig[,2],sep='_'),sig)

freqs <- read.table('../lab.genobaypass_freq',h=F)
freqs$anc <- rowMeans(freqs[,1:2])
dat <- freqs[,c(29,11,12,29,13,14,29,15,16,29,17,18,29,19,29,20,21,29,22,23,29,24,29,25,26,28,27,28)]
colnames(dat) <- c("BSE3-0","BSE3-6","BSE3-10",
				   "BSE4-0","BSE4-6","BSE4-10",
				   "BSE5-0","BSE5-6","BSE5-10",
				   "BSE6-0","BSE6-6","BSE6-10",
				   "BSE7-0","BSE7-6",
				   "BSE8-0","BSE8-6","BSE8-10",
				   "BSE9-0","BSE9-6","BSE9-10",
				   "BSE10-0","BSE10-6",
				   "BSE11-0","BSE11-6","BSE11-10",
				   "BSE12-0","BSE12-6","BSE12-10")

snpdet <- read.table('../lab.snpdet',h=F)[,1:2]
dat <- cbind("GENE.ID"=paste(snpdet[,1],snpdet[,2],sep='_'),dat)

sigdat <- merge(dat,sig,by="GENE.ID")
#polarize by rising allele
for (i in 1:nrow(sigdat)){
	if (sigdat[i,32]<0){
		sigdat[i,2:29] <- 1- sigdat[i,2:29]
			}
		}
# arsine sqrt transform
sigdat[1:nrow(sigdat),2:29] <- 2*asin(sqrt(sigdat[1:nrow(sigdat),2:29]))
sigdat <- sigdat[,1:29]


time <- c("TIME",0,6,10,0,6,10,0,6,10,0,6,10,0,6,0,6,10,0,6,10,0,6,0,6,10,0,6,10)
names(time) <- colnames(sigdat)
sigdat <- rbind(time,sigdat)
bse3 <- sigdat[,1:4]
bse4 <- sigdat[,c(1,5:7)]
bse5 <- sigdat[,c(1,8:10)]
bse6 <- sigdat[,c(1,11:13)]
bse8 <- sigdat[,c(1,16:18)]
bse9 <- sigdat[,c(1,19:21)]
bse11 <- sigdat[,c(1,24:26)]
bse12 <- sigdat[,c(1,27:29)]


write.table(bse3,"bse3.sig_snps.freqs.asinsqrt.txt",row.names=F,quote=F,sep='\t')
write.table(bse4,"bse4.sig_snps.freqs.asinsqrt.txt",row.names=F,quote=F,sep='\t')
write.table(bse5,"bse5.sig_snps.freqs.asinsqrt.txt",row.names=F,quote=F,sep='\t')
write.table(bse6,"bse6.sig_snps.freqs.asinsqrt.txt",row.names=F,quote=F,sep='\t')
write.table(bse8,"bse8.sig_snps.freqs.asinsqrt.txt",row.names=F,quote=F,sep='\t')
write.table(bse9,"bse9.sig_snps.freqs.asinsqrt.txt",row.names=F,quote=F,sep='\t')
write.table(bse11,"bse11.sig_snps.freqs.asinsqrt.txt",row.names=F,quote=F,sep='\t')
write.table(bse12,"bse12.sig_snps.freqs.asinsqrt.txt",row.names=F,quote=F,sep='\t')

#load data
freqdat <- table2eset("bse11.sig_snps.freqs.asinsqrt.txt")

## determine fuzzification
m1 <- mestimate(freqdat)

# determine optimal number of clusters
Dmin(freqdat,m1,crange=2:10,repeats=5,visu=TRUE) #4
cselection(freqdat,m1,crange=2:20,repeats=5,visu=TRUE) #9

#assign clusters to one line
freqdat_cl <- mfuzz(freqdat,c=9,m=m1)
mfuzz.plot(freqdat,cl=freqdat_cl,mfrow=c(2,2),time.labels=c(0,6,10))
O <- overlap(freqdat_cl)
Ptmp <- overlap.plot(freqdat_cl,over=O,thres=0.05)

kclust <- kmeans2(freqdat,4)

#assign membership to other lines


############
# use mean frequencies for each line

sig <- readRDS('../lab.res.new.cmh05.lrt05.RDS')[,c(1,2,20)]
sig <- cbind("GENE.ID"=paste(sig[,1],sig[,2],sep='_'),sig)

dat <- data.frame("Gen0"=rowMeans(freqs[,1:2]),
				"Gen6"=rowMeans(freqs[,c(11,13,15,17,19,20,22,24,25,27)]),
				"Gen10"=rowMeans(freqs[,c(12,14,16,18,21,23,26,28)]))
snpdet <- read.table('../lab.snpdet',h=F)[,1:2]
dat <- cbind("GENE.ID"=paste(snpdet[,1],snpdet[,2],sep='_'),dat)

sigdat <- merge(dat,sig,by="GENE.ID")

#polarize by rising allele
for (i in 1:nrow(sigdat)){
	if (sigdat[i,7]<0){
		sigdat[i,2:4] <- 1- sigdat[i,2:4]
			}
		}

#transform
sigdat[1:nrow(sigdat),2:4] <- 2*asin(sqrt(sigdat[1:nrow(sigdat),2:4]))
sigdat <- sigdat[,1:4]

# calculate divergence from ancestor
sigdat$Gen10 <- sigdat$Gen10 - sigdat$Gen0
sigdat$Gen6 <- sigdat$Gen6 - sigdat$Gen0
sigdat$Gen0 <- sigdat$Gen0 - sigdat$Gen0


time <- c("TIME",0,6,10)
names(time) <- colnames(sigdat)
sigdat <- rbind(time,sigdat)
write.table(sigdat,"mean.sig_snps.div.asinsqrt.txt",row.names=F,quote=F,sep='\t')

#load data
freqdat <- table2eset("mean.sig_snps.div.asinsqrt.txt")

## determine fuzzification
m1 <- mestimate(freqdat)

# determine optimal number of clusters
Dmin(freqdat,m1,crange=2:30,repeats=5,visu=TRUE) #6-7?
cselection(freqdat,m1,crange=2:30,repeats=5,visu=TRUE) #10

freqdat_cl <- mfuzz(freqdat,c=6,m=m1)
mfuzz.plot(freqdat,cl=freqdat_cl,mfrow=c(3,2),time.labels=c(0,6,10))


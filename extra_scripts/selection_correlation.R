# measure correlation in some statistic for different window sizes

library(windowscanr)
library(ggplot2)
library(dplyr)

dat <- readRDS('lab.res.new.RDS')
info <- read.table('lab.snpdet_vs_eaff_position.txt',h=T)
colnames(info) <- c("Pseudo_Transcript", "Pseudo_Position", "Scaffold", "Scaf_Pos")
dat_eaff <- merge(dat,info,by=c("Pseudo_Transcript", "Pseudo_Position"))
dat_eaff_sorted <- dat_eaff[order(dat_eaff$Scaffold,dat_eaff$Scaf_Pos),]

cor_func <- function(mat){
	snpdist <- dist(mat[,1])
	seldist <- dist(mat[,2])
	res <- cor(snpdist,seldist)
	return(res)
	}


winsizes <- c(5,10,50,100)

lrt_matrix <- data.frame("Winsize"=c(),"Distance"=c())
cmh_matrix <- data.frame("Winsize"=c(),"Distance"=c())


for (i in 1:length(winsizes)){
	winsize <- winsizes[i]
	print(winsize)
	winres <- winScan(x = dat_eaff_sorted, 
					 groups = "Pseudo_Transcript", 
					 position = NULL, 
					 values = c("CMH_sel","LRT_pval"), 
					 win_size = winsize,
					 win_step = winsize,
					 funs = c("function (x) mean(dist(-log10(x)))")
					 )
	lrt_res <- data.frame("Winsize"=rep(as.factor(winsize),nrow(winres)),"Distance"=winres[,6])
	cmh_res <- data.frame("Winsize"=rep(as.factor(winsize),nrow(winres)),"Distance"=winres[,8])
	lrt_matrix <- rbind(lrt_matrix,lrt_res)
	cmh_matrix <- rbind(cmh_matrix,cmh_res)

	}



### define window by position, calculate correlation between snp distance and selection distance
## physical distance by lrt distance

randsnps <- sample_n(dat_eaff_sorted,1156)
cmh05 <- readRDS('lab.res.new.cmh05.RDS')
sig <- filter(cmh05,LRT_qval<0.05)
sig_eaff <- merge(sig,info,by=c("Pseudo_Transcript", "Pseudo_Position"))

winsize <- 2000000

sigdat <- data.frame()
for (i in 1:nrow(sig)){
	print(i)
	rowdat <- sig_eaff[i,]
	pos <- rowdat$Scaf_Pos
	scaf <- rowdat$Scaffold
	snps <- filter(dat_eaff_sorted,Scaffold==scaf & Scaf_Pos>(pos-winsize/2) & Scaf_Pos<(pos+winsize/2))
	sel <- -log10(filter(snps,Scaf_Pos==pos)$CMH_sel)
	#sel <- filter(snps,Pseudo_Position==pos)$LRT_pval
	snpdist <- abs(pos-snps$Scaf_Pos)
	snpsel <- -log10(snps$CMH_sel)
	#snpsel <- snps$LRT_pval
	#seldist <- abs(2*asin(sqrt(sel))-2*asin(sqrt(snpsel)))
	#seldist <- abs(sel-snpsel)
	res <- cbind(snpdist,snpsel)
	sigdat <- rbind(sigdat,res)
}

randat <- data.frame()
for (i in 1:nrow(randsnps)){
	rowdat <- randsnps[i,]
	pos <- rowdat$Pseudo_Position
	scaf <- rowdat$Pseudo_Transcript
	snps <- filter(dat_eaff_sorted,Pseudo_Transcript==scaf & Pseudo_Position>(pos-winsize/2) & Pseudo_Position<(pos+winsize/2))
	sel <- -log10(filter(snps,Pseudo_Position==pos)$LRT_pval)
	#sel <- filter(snps,Pseudo_Position==pos)$LRT_pval
	snpdist <- abs(pos-snps$Pseudo_Position)
	snpsel <- -log10(snps$LRT_pval)
	#snpsel <- snps$LRT_pval
	#seldist <- abs(2*asin(sqrt(sel))-2*asin(sqrt(snpsel)))
	seldist <- abs(sel-snpsel)
	res <- cbind(snpdist,seldist)
	randat <- rbind(randat,res)
}

sigdat <- cbind(sigdat,Group=rep("Selected",nrow(sigdat)))
randat <- cbind(randat,Group=rep("Non-Selected",nrow(randat)))

plotdat <- rbind(sigdat,randat)
plotdat <- filter(plotdat,snpdist>0)

p <- ggplot(sigdat, aes(x = snpdist, y = seldist) )
p + stat_smooth(method="lm",formula = y ~ log10(x),se = TRUE) + theme_classic() +xlim(c(1,200))

## physical distance by lrt stat

randsnps <- sample_n(dat_eaff_sorted,1156)
cmh05 <- readRDS('lab.res.new.cmh05.RDS')
#sig <- filter(cmh05,LRT_qval<0.05)
sig <- readRDS('lab.res.new.cmh05.lrt05.RDS')
sig_eaff <- merge(sig,info,by=c("Pseudo_Transcript", "Pseudo_Position"))

winsize <- 2000000

sigdat <- data.frame()
for (i in 1:nrow(sig)){
	print(i)
	rowdat <- sig_eaff[i,]
	pos <- rowdat$Pseudo_Position
	scaf <- rowdat$Pseudo_Transcript
	snps <- filter(dat_eaff_sorted,Pseudo_Transcript==scaf & Pseudo_Position>(pos-winsize/2) & Pseudo_Position<(pos+winsize/2))
	sel <- -log10(filter(snps,Pseudo_Position==pos)$LRT_pval)
	#sel <- abs(filter(snps,Pseudo_Position==pos)$estS.lm)
	snpdist <- abs(pos-snps$Pseudo_Position)
	snpsel <- -log10(snps$LRT_pval)
	#snpsel <- abs(snps$estS.lm)
	res <- cbind(snpdist,snpsel)
	sigdat <- rbind(sigdat,res)
}

sigdat <- filter(sigdat,snpdist>0)
p <- ggplot(sigdat, aes(x = log10(snpdist), y = snpsel))
p +stat_smooth(method="lm",formula=y~x,se = TRUE) + theme_classic() +xlim(c(0,100000))
p + geom_point()

## combined plot
sigdatDist <- cbind(sigdatDist,Group=rep("Distance",nrow(sigdatDist)))
colnames(sigdatDist)[2] <- "snpsel"
sigdat <- cbind(sigdat,Group=rep("Pval",nrow(sigdat)))

plotdat <- rbind(sigdatDist,sigdat)
plotdat <- filter(plotdat,snpdist>0)

p <- ggplot(plotdat, aes(x = snpdist, y = snpsel, color=Group))
p +stat_smooth(method="lm",formula=y~log(x),se = TRUE) + theme_classic() +xlim(c(0,100))


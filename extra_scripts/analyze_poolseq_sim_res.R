library(dplyr)
library(qvalue)

fold_AF <- function(AF){
	if (AF > 0.5 ){
		folded_AF <- 1-AF
		return(folded_AF)
		} else{
		return(AF)
		}	
	}

lrt_res <- readRDS('lrt_res.RDS')
cmh_res <- readRDS('cmh_res.RDS')
freqdat <- read.table('unifsel.10ksel_90knon.freq',h=T)
AF <- rowMeans(freqdat[,c(1,4,7,10,13,16,19,22,25,28,31,34,37,40)])
MAF <- sapply(AF,fold_AF)

res <- data.frame("AF"=AF,"MAF"=MAF,"LRT_pval"=lrt_res,"CMH_pval"=cmh_res[,2], "SNP"=1:10000)
saveRDS(res,'sim.sel_results.RDS')


res$CMH_qval <- qvalue(res$CMH_pval)$qvalues
cmh_sig <- filter(res,CMH_qval<0.05)
cmh_sig$LRT_qval <- qvalue(cmh_sig$LRT_pval)$qvalues
sig <- filter(cmh_sig, LRT_qval<0.05)

nrow(sig)
#TPR
nrow(filter(sig,SNP<=1000)) / 1000
#FPR
nrow(filter(sig,SNP>=1001)) / 9000

mean(sig$MAF)
mean(res[1001:10000,2])

cor.test(sig$MAF,sig$LRT_pval, method="spearman")
cor.test(sig$MAF,sig$CMH_pval, method="spearman")

##########
### analyze AF change vs power
##########

dat <- readRDS('full_data_lmm.raw.sim.10ksel_90knon.RDS')
lrt_res <- readRDS('lrt_res.RDS')
cmh_res <- readRDS('cmh_res.RDS')[,2]

res_all <- as.data.frame(cbind(lrt_res,cmh_res,SNP=1:100000))
res_all$cmh_qval <- qvalue(res_all$cmh_res)$qvalues
res_all$lrt_qval <- qvalue(res_all$lrt_res)$qvalues
res_cmh <- filter(res_all, cmh_qval<0.05)
res_cmh$lrt_qval <- qvalue(res_cmh$lrt_res)$qvalues
sig <- filter(res_cmh,lrt_qval<0.05)
sig_lrt_only <- filter(res_all,lrt_qval<0.05)

af_change6 <- c()
af_change10 <- c()

for (i in 1:10000){
	print(i)
	dt <- filter(dat,SNP==i)
	gen6 <- filter(dt, Treat=="Treatment" & Generation==6)$value
	gen10 <- filter(dt, Treat=="Treatment" & Generation==10)$value
	##gen6
	af_change6 <- c(af_change6,abs(mean(gen6)))
	af_change10 <- c(af_change10,abs(mean(gen10)))
	}

res <- data.frame("SNP"=1:10000,"LRT_pval"=lrt_res[1:10000],"CMH_pval"=cmh_res[1:10000],"Gen6"=af_change6,"Gen10"=af_change10)

af05 <- filter(res,Gen10 <=0.05)
nrow(merge(af05,sig,by="SNP")) / nrow(af05)

af10 <- filter(res,Gen10 > 0.05 & Gen10 <=0.1)
nrow(merge(af10,sig,by="SNP")) / nrow(af10)

af15 <- filter(res,Gen10 > 0.1 & Gen10 <=0.15)
nrow(merge(af15,sig,by="SNP")) / nrow(af15)

af20 <- filter(res,Gen10 > 0.15 & Gen10 <=0.20)
nrow(merge(af20,sig,by="SNP")) / nrow(af20)



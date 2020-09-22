library(qvalue)

fold_AF <- function(AF){
	ifelse (AF > 0.5, 1-AF, AF)
		}

######## real data results

#read in snp details
snpdet <- read.table('lab.snpdet',h=F)

# lrt results
lrt_res <- readRDS('trtxgen_vs_gen.ranBeaker.neff.lrt.pvals.RDS')
qobj <- qvalue(lrt_res)
lrt_qval <- qobj$qvalues

# cmh results
cmh_res <- readRDS('ACER_LLS_runs/acer.cmh.Ne1000.sel.RDS')[,2]

# t test results
#lm_res <- readRDS('basic_lm.trans.neff.pvals.RDS')

#levene results
#lv_res <- readRDS('basic_lv.trans.neff.pvals.RDS')

## put together data frame
res <- data.frame("Pseudo_Transcript"=snpdet[,1],"Pseudo_Position"=snpdet[,2],"LRT_pval"=lrt_res,"CMH_pval"=cmh_res)

######### wild population baypass results
snpdet_wild <- read.table('wild/wild.noPBE.snpdet',h=F)
fst <- readRDS('wild/wild.noPBE.fst.RDS')$snp.FST
fst <- fst[1:655358]
XtX <- readRDS('wild/XtX_res.means.RDS')
XtX_6pops <- readRDS('wild/XtX.6pops.mean.RDS')
std_eBP_6pops <- readRDS('wild/std_eBP.6pops.medians.RDS')$std_eBP
std_eBP <- readRDS('wild/std_eBP.medians.RDS')$std_eBP
std_eBP_bin <- readRDS('wild/std_eBP.bin.medians.RDS')$std_eBP
#pearson_cor <- read.table('wild/wild_std1_summary_betai_reg.out',h=T)$M_Pearson
wild_AF <- read.table('wild/wild.noPBE.genobaypass_freq',h=F)
mean_AF <- rowMeans(wild_AF)

mean_saline_AF <- rowMeans(wild_AF[,c(1,2,3,6,8)]) #above 5psu
mean_fresh_AF <- rowMeans(wild_AF[,c(4,5,7,9,10,11,12)]) #below 5psu
mean_North_Sea_AF <- rowMeans(wild_AF[,c(6,8,11)])
mean_North_Baltic_AF <- rowMeans(wild_AF[,c(1,2,3,4,5,9,10,12)])
mean_Kiel_AF <- wild_AF[,7]

wild_res <- data.frame("Transcript"=snpdet_wild[,1],"Position"=snpdet_wild[,2],"Fst"=fst,"XtX"=XtX$XtX,"std_eBP"=std_eBP,mean_AF,mean_saline_AF,mean_fresh_AF,
				std_eBP_bin,mean_North_Sea_AF,mean_North_Baltic_AF,mean_Kiel_AF,"XtX_6pops"=XtX_6pops$XtX, "std_eBP_6pops"=std_eBP_6pops)

#### get wild stats for lab snps

res_all <- merge(res,wild_res,by=c("Pseudo_Transcript","Pseudo_Position"),all.x=T)
saveRDS(res_all,'lab.res.RDS')
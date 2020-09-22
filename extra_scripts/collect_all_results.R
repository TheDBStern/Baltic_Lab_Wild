library(qvalue)

fold_AF <- function(AF){
	ifelse (AF > 0.5, 1-AF, AF)
		}

######## read in CMH and LRT results

#read in snp details
snpdet <- read.table('lab.snpdet',h=F)

# lrt results
lrt_res <- readRDS('trtxgen_vs_gen.ranBeaker.neff.lrt.pvals.RDS')

# cmh results
cmh_res <- readRDS('ACER_LLS_runs/acer.cmh.Ne1000.sel.RDS')[,2]
cmh_qobj <- qvalue(cmh_res)
cmh_qval <- cmh_qobj$qvalues

# selection coef
selCoef <- readRDS('estS.lm.logit.RDS')

## starting AFS
geno <- read.table("lab.genobaypass_freq",h=F)[,1:2]
T0_AF <- rowMeans(geno)
T0_AF_folded <- fold_AF(T0_AF)

## put together data frame
res <- data.frame("Pseudo_Transcript"=snpdet[,1],"Pseudo_Position"=snpdet[,2],"T0_AF"=T0_AF,"T0_AF_folded"=T0_AF_folded,"LRT_pval"=lrt_res,"CMH_pval"=cmh_res,"CMH_qval"=cmh_qval,"selCoef"=selCoef)
res$T0_AF_rising <- ifelse(res$selCoef<0, 1-res$T0_AF, res$T0_AF)


######### wild results 
snpdet_wild <- read.table('wild/wild.noPBE.snpdet',h=F)
fst <- readRDS('wild/wild.noPBE.fst.RDS')$snp.FST
fst <- fst[1:655358]
XtX <- readRDS('wild/XtX_res.means.RDS')
beta_all <- readRDS('wild/baypass_runs/12pops/std_beta.mean.RDS')[,3]
beta_baltic <- readRDS('wild/baypass_runs/12pops/std_beta.Baltic.mean.RDS')[,3]

#pearson_cor <- read.table('wild/wild_std1_summary_betai_reg.out',h=T)$M_Pearson
wild_AF <- read.table('wild/wild.noPBE.genobaypass_freq',h=F)
mean_AF <- rowMeans(wild_AF)

mean_saline_AF <- rowMeans(wild_AF[,c(1,2,3,6,8)]) #above 5psu
mean_fresh_AF <- rowMeans(wild_AF[,c(4,5,7,9,10,11,12)]) #below 5psu
mean_North_Sea_AF <- rowMeans(wild_AF[,c(6,8,11)])
mean_Baltic_AF <- rowMeans(wild_AF[,c(1,2,3,4,5,9,10,12)])

C2 <- readRDS('wild/baypass_runs/12pops/wild.noPBE.contrasts.RDS')

wild_res <- data.frame("Pseudo_Transcript"=snpdet_wild[,1],"Pseudo_Position"=snpdet_wild[,2],"Fst"=fst,"XtX"=XtX$XtX,
						mean_saline_AF,mean_fresh_AF,mean_North_Sea_AF,mean_Baltic_AF,
						beta_all,beta_baltic,
						"C2_All"=C2$C2_All,"C2_Baltic"=C2$C2_N_Baltic,"C2_N_Sea"=C2$C2_N_Sea)

#### get wild stats for lab snps

res_all <- merge(res,wild_res,by=c("Pseudo_Transcript","Pseudo_Position"),all.x=T)
saveRDS(res_all,'lab.res.new.RDS')
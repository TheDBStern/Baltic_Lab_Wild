## determine number of wild-selected snps among lab-selected SNP
library(dplyr)
library(qvalue)

res_all <- readRDS('lab.res.new.RDS')
res_all$XtXst_qval <- qvalue(1-pchisq(res_all$XtXSt,12))$qvalues
res_all$C2_All_qval <- qvalue(1-pchisq(res_all$C2_All,1))$qvalues
res_all$C2_noKiel_qval <- qvalue(1-pchisq(res_all$C2_noKiel,1))$qvalues
res_all$C2_N_Baltic_qval <- qvalue(1-pchisq(res_all$C2_N_Baltic,1))$qvalues
res_all$C2_N_Sea_qval <- qvalue(1-pchisq(res_all$C2_N_Sea,1))$qvalues


res_cmh05 <- filter(res_all,CMH_sel_qval<0.05)
res_cmh05$LRT_qval <- qvalue(res_cmh05$LRT_pval)$qvalues
sig05 <- filter(res_cmh05,LRT_qval<0.05)
sig01 <- filter(res_cmh05,LRT_qval<0.01)


nrow(filter(sig05,C2_N_Sea > 3.58 & C2_N_Baltic > 3.29 & C2_All > 2.02))
nrow(filter(sig05,C2_N_Sea > 3.58 & C2_N_Baltic > 3.29))
nrow(filter(sig05, C2_N_Baltic > 3.29 & C2_All > 2.02))
nrow(filter(sig05,C2_N_Sea > 3.58& C2_All > 2.02))
nrow(filter(sig05,C2_All > 2.02))
nrow(filter(sig05,C2_N_Sea > 3.58))
nrow(filter(sig05,C2_N_Baltic > 3.29))

wild_res <- read.table('wild/wild.C2_statistics.txt',h=T)
wild_sig <- filter(wild_res,C2_N_Sea > 6.46 | C2_N_Baltic > 5.73 | C2_All > 3.7)

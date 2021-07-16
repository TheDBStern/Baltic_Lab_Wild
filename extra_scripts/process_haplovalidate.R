# for stats need size, nsnps, start, stop
# also - median frequency and coverage across all SNPs for all lines

library(dplyr)

hapval_snps <- readRDS('cmh_qval05.haplotype_blocks.rds')$dominant_haplotypes
info <- read.table('../../lab.snpdet_vs_dovetail_red_positions.txt',h=T)
freq <- read.table('../../lab.snpdet_genobaypass_freq',h=F)
colnames(freq)[1:2] <- c("Pseudo_Transcript","Pseudo_Position")
cov <- read.table('../../lab.snpdet_genobaypass_cov',h=F)
colnames(cov)[1:2] <- c("Pseudo_Transcript","Pseudo_Position")

# convert freq and cov file to dovetail scaffolds and positions
freq_dovetail <- merge(info,freq, by=c("Pseudo_Transcript","Pseudo_Position"))[,3:32]
colnames(freq_dovetail)[1:2] <- c("chr","pos")
cov_dovetail <- merge(info,cov, by=c("Pseudo_Transcript","Pseudo_Position"))[,3:32]
colnames(cov_dovetail)[1:2] <- c("chr","pos")

#collect stats
hap_stats <- hapval_snps %>% group_by(tag) %>%
            summarize(start=min(pos),stop=max(pos),size=max(pos)-min(pos),nsnps=length(pos))

#generate frequencies and coverage for each haplotype
hapval_snps_freq <- merge(hapval_snps,freq_dovetail, by=c("chr","pos"))
hapval_snps_cov <- merge(hapval_snps,cov_dovetail, by=c("chr","pos"))


med_freq <- hapval_snps_freq %>% group_by(tag) %>%
            summarize(V1=median(V3),
                      V2=median(V4),
                      V3=median(V5),
                      V4=median(V6),
                      V5=median(V7),
                      V6=median(V8),
                      V7=median(V9),
                      V8=median(V10),
                      V9=median(V11),
                      V10=median(V12),
                      V11=median(V13),
                      V12=median(V14),
                      V13=median(V15),
                      V14=median(V16),
                      V15=median(V17),
                      V16=median(V18),
                      V17=median(V19),
                      V18=median(V20),
                      V19=median(V21),
                      V20=median(V22),
                      V21=median(V23),
                      V22=median(V24),
                      V23=median(V25),
                      V24=median(V26),
                      V25=median(V27),
                      V26=median(V28),
                      V27=median(V29),
                      V28=median(V30))

med_cov <- hapval_snps_cov %>% group_by(tag) %>%
                    summarize(V1=median(V3),
                              V2=median(V4),
                                            V3=median(V5),
                                            V4=median(V6),
                                            V5=median(V7),
                                            V6=median(V8),
                                            V7=median(V9),
                                            V8=median(V10),
                                            V9=median(V11),
                                            V10=median(V12),
                                            V11=median(V13),
                                            V12=median(V14),
                                            V13=median(V15),
                                            V14=median(V16),
                                            V15=median(V17),
                                            V16=median(V18),
                                            V17=median(V19),
                                            V18=median(V20),
                                            V19=median(V21),
                                            V20=median(V22),
                                            V21=median(V23),
                                            V22=median(V24),
                                            V23=median(V25),
                                            V24=median(V26),
                                            V25=median(V27),
                                            V26=median(V28),
                                            V27=median(V29),
                                            V28=median(V30))

write.table(med_freq[,2:29],'cmh_q05.hapval_freq',quote=F,row.names=F,col.names=F,sep=' ')
write.table(med_cov[,2:29],'cmh_q05.hapval_cov',quote=F,row.names=F,col.names=F,sep=' ')

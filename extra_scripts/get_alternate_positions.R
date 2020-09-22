library(dplyr)

#get pseudoref eaff positions in transcriptome and eaff genome for sig SNPs
# get sig snps
res <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.cmh05.RDS')
sig <- filter(res, LRT_qval<0.05)
sig <- sig[c(1,2,7,18)]
slope <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/sig_snps.rawAF_correlation_coef.RDS')
sig <- cbind(sig,abs(slope))

colnames(sig) <- c("Pseudo_Transcript","Pseudo_Position","LRT_qval","CMH_sel_qval","slope")


#considering also sig in wild
#sig_xtx <- read.table('sig_snps.cmh05.lrt05.XtX95.txt',h=T)[,c(1,2,8,19)]
#sig_ebp <- read.table('sig_snps.cmh05.lrt05.eBP1.3.txt',h=T)[,c(1,2,8,19)]
#sig <- union(sig_xtx,sig_ebp)
#colnames(sig) <- c("Pseudo_Transcript","Pseudo_Position","XtX","eBP")

info_eaff <- read.table('~/Desktop/Baltic_sea_project/pseudoref/varscan/lab.snpdet_vs_eaff_position.txt',h=T)
colnames(info_eaff) <- c("Pseudo_Transcript","Pseudo_Position","Scaffold","Scaf_Position")
info_trans <- read.table('~/Desktop/Baltic_sea_project/pseudoref/varscan/lab.snpdet_vs_transcriptome_positions.txt',h=T)
colnames(info_trans) <- c("Pseudo_Transcript","Pseudo_Position","Tran_Transcript","Tran_Position")

sig_eaff <- merge(sig, info_eaff, by=c("Pseudo_Transcript","Pseudo_Position"))
sig_trans <- merge(sig, info_trans, by=c("Pseudo_Transcript","Pseudo_Position"))

both <- merge(sig_trans, info_eaff, by=c("Pseudo_Transcript","Pseudo_Position"))

## get transcriptome annotation
trans_annotation <- read.delim('~/Desktop/Baltic_sea_project/pseudoref/reference/transcriptome.trinotate_annotation_report.invert.xls',h=T,sep='\t')[,1:5]
colnames(trans_annotation) <- c('Tran_Transcript','All_uniprot_BLASTX','Invert_uniprot_BLASTP','EAFF_i5k_proteins_BLASTP','RefSeq_Invert_BLASTP')

sig_trans_annotation <- merge(both,trans_annotation,all.x=T)

up_info <- read.table('~/Desktop/Baltic_sea_project/pseudoref/varscan/annotation/david_results/sig_snps.cmh05.lrt05.transcriptome_uniprot_genes.invert.txt',sep='\t',h=T)

colnames(up_info) <- c("Invert_uniprot_BLASTP","Uniprot_name","species")
up_info <- up_info[,1:2]

res_up <- merge(sig_trans_annotation,up_info,by="Invert_uniprot_BLASTP",all.x=T)

refseq_info <- read.delim('~/Desktop/Baltic_sea_project/pseudoref/varscan/annotation/sig_snps.cmh05.lrt05.transcriptome_refseq_gene_names.invert.txt',sep='\t',h=T)
colnames(refseq_info) <- c("RefSeq_Invert_BLASTP","RefSeq_name")

res_refseq <- merge(refseq_info,res,by=c("RefSeq_Invert_BLASTP"),all.y=T)

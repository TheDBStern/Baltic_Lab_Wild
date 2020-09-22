library(dplyr)
library(qqman)

res <- readRDS('lab.res.new.RDS')

dovetail_positions <- read.table('lab.snpdet_vs_dovetail_red_positions.txt',h=T)

# determine the main 4 scaffolds
main_scaffolds <- head(dovetail_positions %>% 
  group_by(Dovetail_Scaffold) %>%
  summarise(n_rows = length(Dovetail_Scaffold)) %>%
  arrange(desc(n_rows)),4)$Dovetail_Scaffold

dovetail_positions <- filter(dovetail_positions,Dovetail_Scaffold %in% main_scaffolds)

dovetail_positions$Dovetail_Scaffold <- ifelse(dovetail_positions$Dovetail_Scaffold == main_scaffolds[1],1,dovetail_positions$Dovetail_Scaffold)
dovetail_positions$Dovetail_Scaffold <- ifelse(dovetail_positions$Dovetail_Scaffold == main_scaffolds[2],2,dovetail_positions$Dovetail_Scaffold)
dovetail_positions$Dovetail_Scaffold <- ifelse(dovetail_positions$Dovetail_Scaffold == main_scaffolds[3],3,dovetail_positions$Dovetail_Scaffold)
dovetail_positions$Dovetail_Scaffold <- ifelse(dovetail_positions$Dovetail_Scaffold == main_scaffolds[4],4,dovetail_positions$Dovetail_Scaffold)


res_dovetail <- merge(res,dovetail_positions,by=c("Pseudo_Transcript","Pseudo_Position"))

res_dovetail$snp <- paste(res_dovetail$Dovetail_Scaffold,res_dovetail$Dovetail_Position,sep='_')

res_dovetail$Dovetail_Scaffold <- as.numeric(res_dovetail$Dovetail_Scaffold )

manplot <- manhattan(res_dovetail, chr="Dovetail_Scaffold", bp="Dovetail_Position", snp="snp", p="LRT_pval",
			suggestiveline = -log10(5e-05), genomewideline = -log10(1e-05) )

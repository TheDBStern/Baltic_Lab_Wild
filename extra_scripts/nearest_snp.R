######
## Overlap on the transcript level
######
library(dplyr)

lab_all <- read.table('lab.snpdet',h=F)
res <- readRDS('lab.res.new.cmh05.RDS')
sig <- filter(res,LRT_qval<0.05)

wild_all <- read.table('wild/wild.noPBE.snpdet',h=F)
wild_res <- read.table('wild/wild.C2_statistics.txt',h=T)
wild_sig <- filter(wild_res,C2_N_Sea>6.46)
#wild_sig <- filter(wild_res,C2_All>3.7)
#wild_sig <- filter(wild_res,C2_N_Baltic>5.73)

lab_trans <- unique(sig[,1])
wild_trans <- unique(wild_sig[,1])

length(intersect(lab_trans,wild_trans)) / length(unique(sig[,1]))

overlap <- c()
lab_total <- c()
wild_total <- c()

for (i in 1:1000) {
	lab <- sample_n(lab_all,nrow(sig))
	wild <- sample_n(wild_all,nrow(wild_sig))
	#lab_total <- c(lab_total,length(unique(lab[,1])))
	#wild_total <- c(wild_total,length(unique(wild[,1])))
	res <- length(intersect(unique(lab[,1]),unique(wild[,1]))) / length(unique(lab[,1]))
	overlap <- c(overlap,res)
	}
	
######
## Median distance to nearest SNP in other dataset
######
info <- read.table('lab.snpdet_vs_eaff_position.txt',h=T)
colnames(info)[1:2] <- c("Pseudo_Transcript", "Pseudo_Position")
sig_eaff <- merge(sig,info,by=c("Pseudo_Transcript", "Pseudo_Position"))
sig_eaff <- sig_eaff[,20:21]

wild_sig <- wild_sig[,11:12]
colnames(wild_sig) <- colnames(sig_eaff)

lab_all <- read.table('lab.snpdet_eaff_positions',h=F)
wild_all <- read.table('wild/wild.noPBE.snpdet_trascriptome_eaff_positions',h=T)
wild_all <- wild_all[,5:6]

nearest_dist <- function(row,comp_list){
	scaf <- as.character(row[1])
	pos <- as.numeric(row[2])
	scaf_dat <- comp_list %>% filter(comp_list[,1]==scaf)
	dists <- pos - scaf_dat[,2]
	return(min(abs(dists)))
	}

lab_nearest_wild <- apply(sig_eaff,1, function(x) nearest_dist(x,wild_sig))
lab_nearest_wild <- lab_nearest_wild[!is.infinite(lab_nearest_wild)]

ran_meds <- c()
for (i in 1:1000) {
	print(i)
	lab <- sample_n(lab_all,nrow(sig))
	wild <- sample_n(wild_all,nrow(wild_sig))
	res <- apply(lab,1, function(x) nearest_dist(x,wild))
	res <- res[!is.infinite(res)]
	ran_meds <- c(ran_meds,median(res))
	print(c(mean(res),median(res)))
	}

# median(lab_nearest_wild)
  284
## permute lab, fix wild
ran_meds_shuf_lab <- ran_meds
   2.5%    97.5% 
295.0000 953.0875
## permute wild, fix lab
ran_meds_shuf_wild <- ran_meds
   2.5%   97.5% 
 89.975 144.000 

## permute both
ran_meds_shuf_both <- ran_meds
 2.5% 97.5% 
   58    85 

##########
## Degree of clustering in for sig SNPs
##########

nearest_dist_self <- function(row,comp_list){
	scaf <- as.character(row[1])
	pos <- as.numeric(row[2])
	scaf_dat <- comp_list %>% filter(comp_list[,1]==scaf)
	dists <- pos - scaf_dat[,2]
	dists <- dists[which(dists!=0)]
	return(min(abs(dists)))
	}

samp <- sample_n(wild_all,nrow(wild_sig))

wild_nearest <- apply(samp,1, function(x) nearest_dist_self(x,samp))



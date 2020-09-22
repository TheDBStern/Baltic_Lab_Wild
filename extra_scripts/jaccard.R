## calculate jaccard index 
library(ggplot2)
library(dplyr)
dat <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/full_data_lmm.rawAF.RDS')
snpdet <- read.table('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.snpdet')
dat$Pseudo_Transcript <- rep(snpdet[,1],26)
dat$Pseudo_Position <- rep(snpdet[,2],26)

#res_all <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.new.RDS')
sig <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.new.cmh05.lrt05.RDS')
#non_sig <- filter(res_all, LRT_pval > max(sig$LRT_pval) | CMH_sel_qval >minafc.05)
#non_sig <- sample_n(non_sig,10000)

sig_dat <- merge(dat,sig,by=c("Pseudo_Transcript","Pseudo_Position"))

# for SNPs on same transcript, get mean frequency and selection coef in each line
sig_dat <- sig_dat %>%
	group_by(Pseudo_Transcript,Beaker,Generation) %>%
		summarize(value = mean(value),SNP=first(SNP),Treat=first(Treat),estS_lm=mean(estS_lm))

BSE3_6 <- c()
BSE3_10 <- c()
BSE4_6 <- c()
BSE4_10 <- c()
BSE5_6 <- c()
BSE5_10 <- c()
BSE6_6 <- c()
BSE6_10 <- c()
BSE7_6 <- c()
BSE8_6 <- c()
BSE8_10 <- c()
BSE9_6 <- c()
BSE9_10 <- c()
BSE10_6 <- c()
BSE11_6 <- c()
BSE11_10 <- c()
BSE12_6 <- c()
BSE12_10 <- c()


for (i in 1:length(unique(sig_dat$SNP))){
	snpnum <- unique(sig_dat$SNP)[i]
	print(snpnum)
	dt <- filter(sig_dat,SNP==snpnum)
	#dt <- sig_dat[which(sig_dat$SNP==snpnum)]
	gen6 <- filter(dt, Treat=="Treatment" & Generation==6)$value
	gen10 <- filter(dt, Treat=="Treatment" & Generation==10)$value
	meanafc <- dt$estS_lm[[1]]
	minafc <- 0.05
	##gen6
	if (meanafc*gen6[1]>0 & abs(gen6[1])>minafc){
		BSE3_6 <- c(BSE3_6,dt$SNP[1])
		}
	if (meanafc*gen6[2]>0 & abs(gen6[2])>minafc){
		BSE4_6 <- c(BSE4_6,dt$SNP[1])
		}
	if (meanafc*gen6[3]>0 & abs(gen6[3])>minafc){
		BSE5_6 <- c(BSE5_6,dt$SNP[1])
		}
	if (meanafc*gen6[4]>0 & abs(gen6[4])>minafc){
		BSE6_6 <- c(BSE6_6,dt$SNP[1])
		}
	if (meanafc*gen6[5]>0 & abs(gen6[5])>minafc){
		BSE7_6 <- c(BSE7_6,dt$SNP[1])
		}
	if (meanafc*gen6[6]>0 & abs(gen6[6])>minafc){
		BSE8_6 <- c(BSE8_6,dt$SNP[1])
		}
	if (meanafc*gen6[7]>0 & abs(gen6[7])>minafc){
		BSE9_6 <- c(BSE9_6,dt$SNP[1])
		}
	if (meanafc*gen6[8]>0 & abs(gen6[8])>minafc){
		BSE10_6 <- c(BSE10_6,dt$SNP[1])
		}
	if (meanafc*gen6[9]>0 & abs(gen6[9])>minafc){
		BSE11_6 <- c(BSE11_6,dt$SNP[1])
		}
	if (meanafc*gen6[10]>0 & abs(gen6[10])>minafc){
		BSE12_6 <- c(BSE12_6,dt$SNP[1])
		}
	##gen10
	if (meanafc*gen10[1]>0 & abs(gen10[1])>minafc){
		BSE3_10 <- c(BSE3_10,dt$SNP[1])
		}
	if (meanafc*gen10[2]>0 & abs(gen10[2])>minafc){
		BSE4_10 <- c(BSE4_10,dt$SNP[1])
		}
	if (meanafc*gen10[3]>0 & abs(gen10[3])>minafc){
		BSE5_10 <- c(BSE5_10,dt$SNP[1])
		}
	if (meanafc*gen10[4]>0 & abs(gen10[4])>minafc){
		BSE6_10 <- c(BSE6_10,dt$SNP[1])
		}
	if (meanafc*gen10[5]>0 & abs(gen10[5])>minafc){
		BSE8_10 <- c(BSE8_10,dt$SNP[1])
		}
	if (meanafc*gen10[6]>0 & abs(gen10[6])>minafc){
		BSE9_10 <- c(BSE9_10,dt$SNP[1])
		}
	if (meanafc*gen10[7]>0 & abs(gen10[7])>minafc){
		BSE11_10 <- c(BSE11_10,dt$SNP[1])
		}
	if (meanafc*gen10[8]>0 & abs(gen10[8])>minafc){
		BSE12_10 <- c(BSE12_10,dt$SNP[1])
		}
	}

gen6_pops <- list(BSE3_6,BSE4_6,BSE5_6,BSE6_6,BSE7_6,BSE8_6,BSE9_6,BSE10_6,BSE11_6,BSE12_6)
gen10_pops <- list(BSE3_10,BSE4_10,BSE5_10,BSE6_10,BSE8_10,BSE9_10,BSE11_10,BSE12_10)

gen6_jac <- c()
gen10_jac <- c()

for (pop1 in 1:10){
	for (pop2 in 1:10){
		if (pop1 != pop2){
			jac <- length(intersect(gen6_pops[[pop1]],gen6_pops[[pop2]])) / length(union(gen6_pops[[pop1]],gen6_pops[[pop2]]))
			gen6_jac <- c(gen6_jac,jac)
			}
		}
	}


for (pop1 in 1:8){
	for (pop2 in 1:8){
		if (pop1 != pop2){
			jac <- length(intersect(gen10_pops[[pop1]],gen10_pops[[pop2]])) / length(union(gen10_pops[[pop1]],gen10_pops[[pop2]]))
			gen10_jac <- c(gen10_jac,jac)
			}
		}
	}

gen6_jac_sig <- gen6_jac
gen10_jac_sig <- gen10_jac


#### for simulated data
simdat <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/sims/testing_parallelism/full_data_lmm.sim.1156sel.raw.rep1.RDS')
sig_simdat <- filter(simdat,SNP<=1156)

## loop through all files, calculate jaccard, put into table
plotdat <- data.frame("Group"=c(),"Gen"=c(),"Jaccard"=c())

files <- list.files(path = ".", pattern = 'raw.rep')
for (i in 1:length(files)){
	print(i)
	file <- files[i]
	simdat <- readRDS(file)
	res <- calc_jaccard(simdat)
	df <- data.frame("Group"=rep("Asim",180),"Gen"=c(rep("Asix",90),rep("Ten",90)),"Jaccard"=c(res[[1]],res[[2]]))
	plotdat <- rbind(plotdat,df)
	}

calc_jaccard <- function(sig_simdat){

BSE3_6 <- c()
BSE3_10 <- c()
BSE4_6 <- c()
BSE4_10 <- c()
BSE5_6 <- c()
BSE5_10 <- c()
BSE6_6 <- c()
BSE6_10 <- c()
BSE7_6 <- c()
BSE7_10 <- c()
BSE8_6 <- c()
BSE8_10 <- c()
BSE9_6 <- c()
BSE9_10 <- c()
BSE10_6 <- c()
BSE10_10 <- c()
BSE11_6 <- c()
BSE11_10 <- c()
BSE12_6 <- c()
BSE12_10 <- c()


for (i in 1:length(unique(sig_simdat$SNP))){
	snpnum <- unique(sig_simdat$SNP)[i]
	#print(snpnum)
	dt <- filter(sig_simdat,SNP==snpnum)
	anc <- filter(dt, Treat=="Treatment" & Generation==0)$value
	gen6 <- filter(dt, Treat=="Treatment" & Generation==6)$value - anc
	gen10 <- filter(dt, Treat=="Treatment" & Generation==10)$value - anc
	dt_treat <- filter(dt,Treat=="Treatment")
	dt_treat[which(dt_treat$value==0),3] <- 0.01
	dt_treat[which(dt_treat$value==1),3] <- 0.99
	neff <- (100*dt_treat$Coverage-1)/ (100+dt_treat$Coverage)
	mod <- lmer(log(value/(1-value))~Generation + (1|Beaker), weights=neff, data=dt_treat)
  	sum <- summary(mod)
  	meanafc <- coef(sum)[2,1]
  	minafc <- 0.01
	##gen6
	if (meanafc*gen6[1]>0 & abs(gen6[1])>minafc){
		BSE3_6 <- c(BSE3_6,dt$SNP[1])
		}
	if (meanafc*gen6[2]>0 & abs(gen6[2])>minafc){
		BSE4_6 <- c(BSE4_6,dt$SNP[1])
		}
	if (meanafc*gen6[3]>0 & abs(gen6[3])>minafc){
		BSE5_6 <- c(BSE5_6,dt$SNP[1])
		}
	if (meanafc*gen6[4]>0 & abs(gen6[4])>minafc){
		BSE6_6 <- c(BSE6_6,dt$SNP[1])
		}
	if (meanafc*gen6[5]>0 & abs(gen6[5])>minafc){
		BSE7_6 <- c(BSE7_6,dt$SNP[1])
		}
	if (meanafc*gen6[6]>0 & abs(gen6[6])>minafc){
		BSE8_6 <- c(BSE8_6,dt$SNP[1])
		}
	if (meanafc*gen6[7]>0 & abs(gen6[7])>minafc){
		BSE9_6 <- c(BSE9_6,dt$SNP[1])
		}
	if (meanafc*gen6[8]>0 & abs(gen6[8])>minafc){
		BSE10_6 <- c(BSE10_6,dt$SNP[1])
		}
	if (meanafc*gen6[9]>0 & abs(gen6[9])>minafc){
		BSE11_6 <- c(BSE11_6,dt$SNP[1])
		}
	if (meanafc*gen6[10]>0 & abs(gen6[10])>minafc){
		BSE12_6 <- c(BSE12_6,dt$SNP[1])
		}
	##gen10
	if (meanafc*gen10[1]>0 & abs(gen10[1]) >minafc){
		BSE3_10 <- c(BSE3_10,dt$SNP[1])
		}
	if (meanafc*gen10[2]>0 & abs(gen10[2]) >minafc){
		BSE4_10 <- c(BSE4_10,dt$SNP[1])
		}
	if (meanafc*gen10[3]>0 & abs(gen10[3]) >minafc){
		BSE5_10 <- c(BSE5_10,dt$SNP[1])
		}
	if (meanafc*gen10[4]>0 & abs(gen10[4]) >minafc){
		BSE6_10 <- c(BSE6_10,dt$SNP[1])
		}
	if (meanafc*gen10[5]>0 & abs(gen10[5]) >minafc){
		BSE7_10 <- c(BSE7_10,dt$SNP[1])
		}
	if (meanafc*gen10[6]>0 & abs(gen10[6]) >minafc){
		BSE8_10 <- c(BSE8_10,dt$SNP[1])
		}
	if (meanafc*gen10[7]>0 & abs(gen10[7]) >minafc){
		BSE9_10 <- c(BSE9_10,dt$SNP[1])
		}
	if (meanafc*gen10[8]>0 & abs(gen10[8]) >minafc){
		BSE10_10 <- c(BSE10_10,dt$SNP[1])
		}
	if (meanafc*gen10[9]>0 & abs(gen10[9]) >minafc){
		BSE11_10 <- c(BSE11_10,dt$SNP[1])
		}
	if (meanafc*gen10[10]>0 & abs(gen10[10]) >minafc){
		BSE12_10 <- c(BSE12_10,dt$SNP[1])
		}
	}

gen6_pops <- list(BSE3_6,BSE4_6,BSE5_6,BSE6_6,BSE7_6,BSE8_6,BSE9_6,BSE10_6,BSE11_6,BSE12_6)
gen10_pops <- list(BSE3_10,BSE4_10,BSE5_10,BSE6_10,BSE7_10,BSE8_10,BSE9_10,BSE10_10,BSE11_10,BSE12_10)

gen6_jac <- c()
gen10_jac <- c()

for (pop1 in 1:10){
	for (pop2 in 1:10){
		if (pop1 != pop2){
			jac <- length(intersect(gen6_pops[[pop1]],gen6_pops[[pop2]])) / length(union(gen6_pops[[pop1]],gen6_pops[[pop2]]))
			gen6_jac <- c(gen6_jac,jac)
			}
		}
	}


for (pop1 in 1:10){
	for (pop2 in 1:10){
		if (pop1 != pop2){
			jac <- length(intersect(gen10_pops[[pop1]],gen10_pops[[pop2]])) / length(union(gen10_pops[[pop1]],gen10_pops[[pop2]]))
			gen10_jac <- c(gen10_jac,jac)
			}
		}
	}

gen6_jac_sim <- gen6_jac
gen10_jac_sim <- gen10_jac

return(list(gen6_jac,gen10_jac))
}

mean(gen6_jac_sig)
mean(gen6_jac_sim)
wilcox.test(gen6_jac_sig,gen6_jac_sim)$p.value

mean(gen10_jac_sig)
mean(gen10_jac_sim)
wilcox.test(gen10_jac_sig,gen10_jac_sim)$p.value

########
#plots

plotdat <- data.frame("Group"=c(rep("ASim",length(unique(gen6_jac_sim))+length(unique(gen10_jac_sim))),rep("Emp",length(unique(gen6_jac_sig))+length(unique(gen10_jac_sig)))),
						"Gen"=c(rep("Six",length(unique(gen6_jac_sim))),rep("Ten",length(unique(gen10_jac_sim))),rep("Six",length(unique(gen6_jac_sig))),rep("Ten",length(unique(gen10_jac_sig)))), 
					"Jaccard" = c(unique(gen6_jac_sim),unique(gen10_jac_sim),unique(gen6_jac_sig),unique(gen10_jac_sig)))

saveRDS(plotdat,"jaccard_plotdat.min0.05.selstart.RDS")
########


###
### for controls
BSE1_6 <- c()
BSE1_10 <- c()
BSE2_6 <- c()
BSE2_10 <- c()


for (i in 1:length(unique(sig_dat$SNP))){
	snpnum <- unique(sig_dat$SNP)[i]
	print(snpnum)
	dt <- filter(sig_dat,SNP==snpnum)
	#dt <- sig_dat[which(sig_dat$SNP==snpnum)]
	gen6 <- filter(dt, Treat=="Treatment" & Generation==6)$value
	gen10 <- filter(dt, Treat=="Treatment" & Generation==10)$value
	meanafc <- dt$estS[[1]]
	gen6_ctrl <- filter(dt, Treat=="Control" & Generation==6)$value
	gen10_ctrl <- filter(dt, Treat=="Control" & Generation==10)$value
	if (meanafc*gen6_ctrl[1]>0){
		BSE1_6 <- c(BSE1_6,dt$SNP[1])
		}
	if (meanafc*gen6_ctrl[2]>0){
		BSE2_6 <- c(BSE2_6,dt$SNP[1])
		}
	if (meanafc*mean(gen10_ctrl[1:2])>0){
		BSE1_10 <- c(BSE1_10,dt$SNP[1])
		}
	if (meanafc*mean(gen10_ctrl[3:4])>0){
		BSE2_10 <- c(BSE2_10,dt$SNP[1])
		}
	}

gen6_pops <- list(BSE1_6,BSE2_6)
gen10_pops <- list(BSE1_10,BSE2_10)

gen6_jac <- c()
gen10_jac <- c()

for (pop1 in 1:2){
	for (pop2 in 1:2){
		if (pop1 != pop2){
			jac <- length(intersect(gen6_pops[[pop1]],gen6_pops[[pop2]])) / length(union(gen6_pops[[pop1]],gen6_pops[[pop2]]))
			gen6_jac <- c(gen6_jac,jac)
			}
		}
	}


for (pop1 in 1:2){
	for (pop2 in 1:2){
		if (pop1 != pop2){
			jac <- length(intersect(gen10_pops[[pop1]],gen10_pops[[pop2]])) / length(union(gen10_pops[[pop1]],gen10_pops[[pop2]]))
			gen10_jac <- c(gen10_jac,jac)
			}
		}
	}

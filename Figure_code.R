########
## Working R code for making figures for
## Stern DB, Diaz JA, and CE Lee. Parallel polygenic adaptation driven by epistasis in laboratory and wild populations of a Baltic Sea Copepod
########

########## Fig 1a
#plot allele frequencies by generation, treatment vs control ############

library(ggplot2)
library(dplyr)
dat <- readRDS('full_data.rawFreqs.RDS') # this file can be generated from prep_lmm.R by removing the arcsin transformations and calculating the divergence from the starting population
bse0_tmp <- filter(dat,Beaker=="BSE-0")
bse0_tmp$Treat <- rep("Treatment",nrow(bse0_tmp))
bse0b_tmp <- filter(dat,Beaker=="BSE-0B")
bse0b_tmp$Treat <- rep("Treatment",nrow(bse0b_tmp))
dat <- rbind(dat,bse0_tmp)
dat <- rbind(dat,bse0b_tmp)
dat$folded <- sapply(dat$value,fold_AF)


res <- readRDS('data/lab.all.RDS')
sig <- readRDS('data/lab.sig.RDS')
colnames(sig)[1:2] <- c('Transcript','Position')
sig_dat <- merge(dat,sig,by=c('Transcript','Position'))

## polarize so show SNP rising in freq in selected lines
sig_dat <- as.data.frame(sig_dat)

snpnums <- unique(sig_dat$SNP)

polarized <- data.frame()
for (i in 1:length(snpnums)){
	snpnum <- snpnums[i]
	print(snpnum)
	dt <- filter(sig_dat,SNP==snpnum)
	slope <- dt$selCoef[[1]]
	if (slope < 0 ){
		dt$value <- (1- dt$value)
		}
	polarized <- rbind(polarized, dt)
	}

#plot polarized AF
p <- ggplot(polarized, aes(x = Generation, y = value, color=Treat,fill=Treat))
p +   stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=1, position=position_dodge(width=1)) +
  stat_summary(fun=mean, geom="point", position=position_dodge(width=1)) +
	stat_smooth(method="lm",fullrange = FALSE) +
	scale_color_manual(values=c("grey","#DDAFFA")) +
	scale_fill_manual(values=c("grey","#DDAFFA")) +
	theme_classic()


########## Fig 1c
#plot polarized AFS for starting pop and gen 10
plotdat <- filter(polarized,Treat=="Treatment")
plotdat <- filter(plotdat,Generation==0 | Generation==10)

p <- ggplot(plotdat, aes(x = value, color=factor(Generation),fill=factor(Generation)))
p + geom_density(adjust=0.2,alpha=0.1) +
	scale_color_manual(values=c("#f8766d","#7501BC")) +
	scale_fill_manual(values=c("#f8766d","#7501BC")) +
	theme_classic()

# get starting and ending freqs for every SNP
gen10 <- plotdat %>%
		filter(Generation==10) %>%
		group_by(Transcript,Position) %>%
		summarize(T10_treat_AF_rising=mean(value))

########### Figure 2a #################
library(dplyr)
library(ggbiplot)
library(tcR)
# PCA for sig SNP
freqs <- read.table("lab.genobaypass_freq", h=F) #this file is generated with baypass2freqs_cov.py
snpdet <- read.table("data/lab.snpdet", h=F)
dat <- cbind(snpdet[,1:2],freqs)
colnames(dat) <- c("Pseudo_Transcript","Pseudo_Position","BSE-0","BSE-0B","BS3C","BS4C","BSE-1-T1","BSE-1-T2",
				"BSE-1-T2B","BSE-2-T1","BSE-2-T2","BSE-2-T2B","BSE-3-T1","BSE-3-T2","BSE-4-T1","BSE-4-T2","BSE-5-T1",
				"BSE-5-T2","BSE-6-T1","BSE-6-T2","BSE-7-T1","BSE-8-T1","BSE-8-T2","BSE-9-T1","BSE-9-T2","BSE-10-T1",
				"BSE-11-T1","BSE-11-T2","BSE-12-T1","BSE-12-T2")

group <- c(rep("A",2),rep("C",8),"T6","T10","T6","T10","T6","T10","T6","T10","T6","T6","T10","T6","T10","T6","T6","T10","T6","T10")
generation <- as.character(c(rep(0,2),rep(20,2),6,10,10,6,10,10,6,10,6,10,6,10,6,10,6,6,10,6,10,6,6,10,6,10))
sig <- readRDS('lab.sig.RDS')

sig_dat <- merge(dat,sig,by=c("Pseudo_Transcript","Pseudo_Position"))[,3:30]
sig_dat.t <- t(sig_dat)


pca_res <- prcomp(sig_dat.t, center = TRUE,scale. = TRUE)

ggbiplot(pca_res,var.axes=FALSE,groups=group,ellipse=TRUE) +
	theme_classic()

euc_res <- pca2euclid(pca_res)

treat_T1_pops <- c(11,13,15,17,19,20,22,24,25,27)
T1 <- euc_res[treat_T1_pops,treat_T1_pops]
mean(T1[lower.tri(T1)])

treat_T2_pops <- c(12,14,16,18,21,23,26,28)
T2 <- euc_res[treat_T2_pops,treat_T2_pops]
mean(T2[lower.tri(T2)])


########### Figure 2b #################
#### Assessing parallelism using pairwise overlap (Jaccard Index)
#
library(ggplot2)
library(dplyr)
library(lme4)
source("extra_scripts/parallelism_functions.R")

# empirical data
# data are formatted to show change from ancestor
dat <- readRDS('full_data_lmm.rawAF.RDS')
snpdet <- read.table('data/lab.snpdet')
dat$Pseudo_Transcript <- rep(snpdet[,1],26)
dat$Pseudo_Position <- rep(snpdet[,2],26)

sig <- readRDS('data/lab.sig.RDS')

sig_dat <- merge(dat,sig,by=c("Pseudo_Transcript","Pseudo_Position"))

####
# jaccard plots under different simulation models
####

#empirical data
jaccard_emp <- calc_jaccard_emp(sig_dat,0.01)

#independent SNP sims (sweeps)
jaccard_ind_6 <- c()
jaccard_ind_10 <- c()

ind_sim_files <- list.files(path = "sims/sweep_sims", pattern = 'sel.raw.rep', full.names=T) #generated from the extra_scripts/simulate_poolseq.R
for (i in 1:length(ind_sim_files)){
	print(i)
	file <- ind_sim_files[i]
	simdat <- readRDS(file)
	res <- calc_jaccard_ind_sims(simdat, 0.01)
	jaccard_ind_6 <- c(jaccard_ind_6,res[[1]])
	jaccard_ind_10 <- c(jaccard_ind_10,res[[2]])
	}

#polygenic shift sims, 30 loci
jaccard_qt30_6 <- c()
jaccard_qt30_10 <- c()

qt30_sim_files <- list.files(path = "sims/qt_sims/",
									pattern = 'loc30.rep', full.names=T)

for (i in 1:length(qt30_sim_files)){
	print(i)
	file <- qt30_sim_files[i]
	simdat <- read.table(file,h=F)
	colnames(simdat) <- c("snp_id","generation","frequency","replicate")
	res <- calc_jaccard_qt_sims(simdat, 0.01)
	jaccard_qt30_6 <- c(jaccard_qt30_6,res[[1]])
	jaccard_qt30_10 <- c(jaccard_qt30_10,res[[2]])
	}

#QT sims, 1156 loci
jaccard_qt1156_6 <- c()
jaccard_qt1156_10 <- c()

qt1156_sim_files <- list.files(path = "sims/qt_sims/",
									pattern = 'loc1156', full.names=T)

for (i in 1:length(qt1156_sim_files)){
	print(i)
	file <- qt1156_sim_files[i]
	simdat <- read.table(file,h=F)
	colnames(simdat) <- c("snp_id","generation","frequency","replicate")
	res <- calc_jaccard_qt_sims(simdat, 0.01)
	jaccard_qt1156_6 <- c(jaccard_qt1156_6,res[[1]])
	jaccard_qt1156_10 <- c(jaccard_qt1156_10,res[[2]])
	}

plotdat <- data.frame("Group"=c(rep("AQT1156",length(jaccard_qt1156_6)),
																	   rep("Bind1156",length(jaccard_ind_6)),
																		 rep("CEmp",length(jaccard_emp[[1]])),
																		 rep("DQT30",length(jaccard_qt30_6)),
																		 rep("AQT1156",length(jaccard_qt1156_10)),
																		 rep("Bind1156",length(jaccard_ind_10)),
																		 rep("CEmp",length(jaccard_emp[[2]])),
																		 rep("DQT30",length(jaccard_qt30_10))),
											"Gen"=c(rep("Asix",length(c(jaccard_qt1156_6,jaccard_ind_6,jaccard_emp[[1]],jaccard_qt30_6))),
															rep("Ten",length(c(jaccard_qt1156_10,jaccard_ind_10,jaccard_emp[[2]],jaccard_qt30_10)))),
											"Jaccard"=c(jaccard_qt1156_6,jaccard_ind_6,jaccard_emp[[1]],jaccard_qt30_6,
																	jaccard_qt1156_10,jaccard_ind_10,jaccard_emp[[2]],jaccard_qt30_10))
plotdat$GroupGen <- paste(plotdat$Group,plotdat$Gen,sep="_")

p <- ggplot(plotdat, aes(x=Gen, y=Jaccard, color=Group))
p + geom_boxplot() +
	geom_point(pch = 21, size=0.5, position = position_jitterdodge(jitter.width =0.1))+
 	scale_color_manual(values=c("#CCCCCC","#ABC1D3","#D206FC","#606060")) + #Light Blue, Orange, Purple,  Dark Blue
 theme_classic()

getstats <- function(x){
					print("Mean")
					print(mean(x))
					print("N")
					print(length(x))
					print("Standard Error")
					print(mean(x) / sqrt(length(x)))
}

getstats(jaccard_qt30_10)

#wilcoxon rank sum tests
wilres <- wilcox.test(Jaccard~Group,
							data=filter(plotdat,
								Gen=="Ten" &
								(Group=="CEmp" | Group=="AQT1156")))
wilres$p.value
wilres$statistic

#increase in jaccard over time
lm_res <- lm(Jaccard~Gen, data=filter(plotdat, Group=="DQT30"))

emp_qt30 <- lm(Jaccard~Gen*Group, data=filter(plotdat, Group=="Bind1156" | Group=="CEmp"))



########### Figure 2c #################
#### Assessing parallelism using the replicate frequency spectrum

# data are formatted to show change from ancestor
dat <- readRDS('full_data_lmm.rawAF.RDS')
snpdet <- read.table('data/lab.snpdet')
dat$Pseudo_Transcript <- rep(snpdet[,1],26)
dat$Pseudo_Position <- rep(snpdet[,2],26)

sig <- readRDS('data/lab.sig.RDS')

sig_dat <- merge(dat,sig,by=c("Pseudo_Transcript","Pseudo_Position"))

## calculate rfs
# creates list where first item is gen 6 and second is gen 10
rfs_emp <- rfs_empirical(sig_dat,0.01)

#independent snp simulations
# Loop through each of the 10 simulations and add to vector
# Data are formatted with raw snp frequencies for each snp, so function estimates s and performed calculation
# Again, the first item of the list is gen 6 and second is gen 10

rfs_ind_sim <- c()

ind_sim_files <- list.files(path = "sims/sweep_sims", pattern = 'raw.rep', full.names=T)
for (i in 1:length(ind_sim_files)){
	print(i)
	file <- ind_sim_files[i]
	simdat <- readRDS(file)
	res <- rfs_ind_sims(simdat, 0.01)
	rfs_ind_sim <- c(rfs_ind_sim,res[[2]])
	}

#quantitative trait (polygenic shift) simulations
# same as above, but perform for 30 qtl and 1156 qtl
rfs_qt30 <- c()

qt30_sim_files <- list.files(path = "sims/qt_sims/",
									pattern = 'loc30', full.names=T)

for (i in 1:length(qt30_sim_files)){
	print(i)
	file <- qt30_sim_files[i]
	simdat <- read.table(file,h=F)
	colnames(simdat) <- c("snp_id","generation","frequency","replicate")
	res <- rfs_qt_sims(simdat, 0.01)
	rfs_qt30 <- c(rfs_qt30,res[[2]])
	}

rfs_qt1156 <- c()

qt1156_sim_files <- list.files(path = "sims/qt_sims/",
									pattern = 'loc1156', full.names=T)

for (i in 1:length(qt1156_sim_files)){
	print(i)
	file <- qt1156_sim_files[i]
	simdat <- read.table(file,h=F)
	colnames(simdat) <- c("snp_id","generation","frequency","replicate")
	res <- rfs_qt_sims(simdat, 0.01)
	rfs_qt1156 <- c(rfs_qt1156,res[[2]])
	}

## rfs plots
plotdat <- data.frame("Group"=c(rep("CEmp",length(rfs_emp[[2]])),rep("BInd_sim",length(rfs_ind_sim)),rep("DQT30",length(rfs_qt30)),rep("AQT1156",length(rfs_qt1156))),
					  "Proportion"=c(rfs_emp[[2]],rfs_ind_sim,rfs_qt30,rfs_qt1156))
p <- ggplot(plotdat, aes(x=Proportion, fill=Group))
p + geom_histogram(aes(y=..density..),binwidth=0.1,position=position_dodge2(reverse = FALSE,preserve ="single"))+
	scale_fill_manual("Group",values=c("#ffd79d","grey","#c98bf7","#eaac59")) +
	theme_classic()

#kolmogorov smirnoff tests
ks.test(filter(plotdat,Group=="CEmp")$Proportion, filter(plotdat,Group=="DQT30")$Proportion)
#emp vs ind SNPs
#D = 0.53452, p-value < 2.2e-16

#emp vs QT1156
#D = 0.98953, p-value < 2.2e-16

#emp vs QT30
#D = 0.12095, p-value = 0.001882


######Figure 3a #######
### Map, Baltic Only ####
require(dplyr); require(RColorBrewer); require(ggplot2)
require(mapdata); require(maptools)
#require(rgdal)

sample_data <- read.csv("wild_sample_data.baltic.csv", header = TRUE)

europe <- map_data("world2")

Euromap <- ggplot() + geom_polygon(data = europe,
                                 aes(x=long, y = lat, group = group),
                                 fill = "#EFEFEF",
                                 color="black") +
  coord_map(projection = "albers", lat0 = 53, lat1 = 66, xlim = c(6,26), ylim = c(53,66))

Euromap +
    geom_point(data=sample_data, aes(x=Longitude, y=Latitude, fill=MAS),
               color = "black", shape=21, size=5.0) +
		scale_fill_gradient(low = "#DDAFFA", high = "#7501BC") +
    theme(panel.background = element_rect(fill = "#EBF5FB"),
    		axis.title.x=element_blank(),
    		axis.title.y=element_blank())

######Figure 3c #######
### MAF distributions ####
#c. Folded AFS for selected vs non-selected SNP

fold_af <- function(af){
		ifelse(af > 0.5, 1-af,af)
}
res_all <- readRDS('data/lab.all.RDS')
sig <- readRDS('data/lab.sig.RDS')
non_sig <- filter(res_all, LRT_pval > max(sig$LRT_pval) | CMH_pval > 0.05)

wild_freqs <- read.table('wild.genobaypass_freq',h=F)
colnames(wild_freqs) <- c("BB1E","BB2E","GBE","HF1E","HF2E","KIE","RG1E","RG2E","STE")
wild_snpdet <- read.table('wild.snpdet',h=F)[,1:2]
colnames(wild_snpdet) <- c("Pseudo_Transcript","Pseudo_Position")
wild_dat <- cbind(wild_snpdet,wild_freqs)

wild_sig <- merge(wild_dat,sig, by=c("Pseudo_Transcript","Pseudo_Position"))
wild_non_sig <- merge(wild_dat,non_sig, by=c("Pseudo_Transcript","Pseudo_Position"))

ksres <- ks.test(fold_af(wild_sig$STE), fold_af(wild_non_sig$STE))
ksres$statistic
ksres$p.value

plotdat <- data.frame("Group"=c(rep("Selected",nrow(wild_sig)),rep("Non_Selected",nrow(wild_non_sig))),
					"MAF"= c(fold_af(wild_sig$BB1E), fold_af(wild_non_sig$BB1E)))


library(plyr); library(dplyr)
library(ggplot2)

cdat <- ddply(plotdat, "Group", summarise, MAF.mean=mean(MAF))
p <- ggplot(plotdat, aes(x = MAF, after_stat(density), color=Group,fill=Group))
p + geom_histogram(binwidth = 0.05, position="dodge")+
	geom_vline(data=cdat, aes(xintercept=MAF.mean,  colour=Group),
               linetype="dashed", size=1) +
    scale_color_manual(values=c("#d6d6d6","#DDAFFA")) +
	scale_fill_manual(values=c("#949494","#c870ff")) +
	theme_classic()

ks.test(filter(plotdat, Group=="Selected")$MAF,filter(plotdat, Group=="Non_Selected")$MAF)

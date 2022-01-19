########
## Working R code for making figures for
## Stern DB, Diaz JA, and CE Lee. Parallel polygenic adaptation driven by epistasis in laboratory and wild populations of a Baltic Sea Copepod
########

########## Fig 1a
#manhattan plot ############
## show LMM and CMH p vals, with sign highlighted in red. Haploblocks underneath
library(RColorBrewer)
library(data.table)
library(qvalue)
library(dplyr)

snpdat <- as.data.table(readRDS('lab.snp_res_dt.RDS'))
cmh <- as.data.table(readRDS('acer.cmh.NeVar.sel.RDS'))
cmh$qvalue <- qvalue(10^(-cmh$score))$qvalues
cmh_sig <- filter(cmh, qvalue<0.05)
cmh_cutoff <- filter(cmh_sig, score==min(cmh_sig$score))$score

lmm <- data.table("chr"=snpdat$Dovetail_Scaffold,
                  "pos"=snpdat$Dovetail_Position,
                  "score"= -log10(snpdat$LRT_pval))
lmm$qvalue <- qvalue(10^(-lmm$score))$qvalues
lmm_sig <- filter(lmm, qvalue<0.05)
lmm_cutoff <- filter(lmm_sig, score==min(lmm_sig$score))$score

blocks <- readRDS('combined.cmh_chisq_lrt.haplotype_blocks.RDS')$dominant_haplotypes

sig <- merge(cmh,blocks,)


##grey boxes for hbs and colored if significant
hb_info <- readRDS("hap_blocks.res.RDS")

#plot
newcol.raw <- unique(c(brewer.pal(8,"Set1"),brewer.pal(7,"Dark2"),brewer.pal(7,"Set2"),brewer.pal(12,"Set3"),brewer.pal(12,"Paired"),brewer.pal(7,"Accent"),brewer.pal(7,"Spectral")))
newcol  <- colorRampPalette(newcol.raw)(max(as.numeric(factor(blocks$tag)),na.rm=TRUE))
## as.numeric removes the non-grouped tags
#col.sub <- newcol[as.numeric(factor(blocks$tag))]
#blocks[,groupcol := col.sub]
j <- "Scz6wRH_1693;HRSCAF=1917"
hb_info.sub <- filter(hb_info, chr==j)
#cmh.sub <- merge(blocks[chr==j,],cmh,by=c("chr","pos"))
#cmh.sub <- filter(cmh,qvalue<0.05)
cmh.sub <- data.frame()
for (hb in 1:nrow(hb_info.sub)){
    dat <- filter(cmh, chr==j & pos >= hb_info.sub[hb,"start"] & pos <= hb_info.sub[hb,"stop"] & qvalue < 0.05)
    dat$groupcol <- newcol[[hb]]
    cmh.sub <- rbind(cmh.sub,dat)
}
#lmm.sub <- merge(blocks[chr==j,],lmm,by=c("chr","pos"))
#lmm.sub <- filter(lmm,qvalue<0.05)
lmm.sub <- data.frame()
for (hb in 1:nrow(hb_info.sub)){
    dat <- filter(lmm, chr==j & pos >= hb_info.sub[hb,"start"] & pos <= hb_info.sub[hb,"stop"] & qvalue < 0.05)
    dat$groupcol <- newcol[[hb]]
    lmm.sub <- rbind(lmm.sub,dat)
}
pdf(paste("test_boxes","_",j,"_haplovalidate.pdf",sep=""), width=8,height=2)
par(mfrow=c(2,1), mar=c(0.1, 4, 0.1, 2))
maxy_cmh <- round(max(as.numeric(cmh[chr==j,score]),na.rm=TRUE),0) + round(max(as.numeric(cmh[chr==j,score]),na.rm=TRUE),0) * 0.02
maxy_lmm <- round(max(as.numeric(lmm[chr==j,score]),na.rm=TRUE),0) + round(max(as.numeric(lmm[chr==j,score]),na.rm=TRUE),0) * 0.02
cmh[chr==j,plot(pos/1000000,as.numeric(score),pch=19,cex=0.2,col="#414547",ylim=c(0,maxy_cmh),ylab="CMH score",cex.axis=1,cex.lab=1,xaxt='n')]
cmh.sub[,points(as.numeric(pos)/1000000,as.numeric(score),cex=0.2,pch=19,col=groupcol)]
for (i in 1:nrow(hb_info.sub)){
    rect(xleft=hb_info.sub[i,"start"]/1000000,ybottom=0,xright=hb_info.sub[i,"stop"]/1000000,ytop=20,
    col=rgb(0.5,0.5,0.5,alpha=0.2),border=NA)
}
abline(h=cmh_cutoff,lty = 2,col="red")
lmm[chr==j,plot(pos/1000000,as.numeric(score),pch=19,cex=0.2,col="#414547",ylim=rev(c(0,maxy_lmm)),xlab="position Mb",ylab="LMM score",cex.axis=1,cex.lab=1,xaxp=c(0,150,5))]
lmm.sub[,points(as.numeric(pos)/1000000,as.numeric(score),cex=0.2,pch=19,col=groupcol)]
for (i in 1:nrow(hb_info.sub)){
    rect(xleft=hb_info.sub[i,"start"]/1000000,ybottom=0,xright=hb_info.sub[i,"stop"]/1000000,ytop=10,
    col=rgb(0.5,0.5,0.5,alpha=0.2),border=NA)
}
abline(h=lmm_cutoff,lty = 2,col="red")
dev.off()

########## Fig 1b
# histogram of allele frequency changes / selection coefficients
library(ggplot2)

#selection coefficients
dat <- readRDS('../data/hap_blocks.res.RDS')

p <- ggplot(dat, aes(selCoef))+
    geom_histogram(color="#7501BC", fill="#7501BC", alpha=0.5) +
    geom_vline(aes(xintercept=mean(selCoef)),
            color="black", linetype="dashed", size=1) +
    theme_classic()
p

#starting frequencies
p <- ggplot(dat, aes(T0_AF))+
    geom_histogram(color="#F0765F", fill="#F0765F", alpha=0.5) +
    geom_vline(aes(xintercept=mean(selCoef)),
            color="black", linetype="dashed", size=1) +
    theme_classic()
p

########## Fig 1c
#plot polarized AFS for starting pop and gen 10
plotdat <- filter(dat,Treat=="Treatment")
plotdat <- filter(dat,Generation==0 | Generation==10)

p <- ggplot(plotdat, aes(x = value, color=factor(Generation),fill=factor(Generation)))
p + geom_density(adjust=0.7,alpha=0.1) +
	scale_color_manual(values=c("#f8766d","#7501BC")) +
	scale_fill_manual(values=c("#f8766d","#7501BC")) +
	theme_classic()

# get starting and ending freqs for every SNP
gen10 <- plotdat %>%
		filter(Generation==10) %>%
		group_by(Transcript,Position) %>%
		summarize(T10_treat_AF_rising=mean(value))

###########
### Allele frequency plot
###########
library(ggplotFL)
library(ggplot2)
library(dplyr)

##### for hbs
dat <- readRDS("hap_blocks.rawFreqs.RDS")

#` make dummy ancestor data for treatment lines
anc <- filter(dat,variable=="Ancestor")
anc$Treat <- "Treatment"
dat <- rbind(dat,anc)

pdat <- dat %>%
        group_by(SNP,Generation,Treat) %>%
        summarise(AF=mean(value))
##simulate data from mean starting frequency
library(poolSeq)
startingAF <- mean(filter(pdat,Generation==0)$AF)
sims <- wf.traj(p0=rep(startingAF,10000),Ne=1750,t=c(0,6,10))
simdat <- data.frame(AF=c(sims[,1],sims[,2],sims[,3]),
                     SNP=rep(1:10000,3),
                    Generation=c(rep(0,10000),rep(6,10000),rep(10,10000)))

p <-  ggplot(filter(pdat,Treat=="Treatment"), aes(y = AF, x = Generation,group=SNP))

p + geom_line(alpha=0.2) +
    geom_point(alpha=0.2) +
    stat_summary(data=filter(pdat,Treat=="Treatment" & (Generation == 0 | Generation == 10)),fun=mean, size=2,colour="#C98BF7", geom="line", aes(group = 1))+
    stat_summary(data=filter(pdat,Treat=="Control" & (Generation == 0 | Generation == 20)),fun=mean, size=2,colour="#F9BE70", geom="line", aes(group = 1))+
    geom_flquantiles(data=simdat,aes(Generation,AF,group=NULL),probs=c(0.01,0.99), fill="blue", alpha=0.1) +
    theme_classic() +
    coord_cartesian(xlim=c(0, 10))

########### Figure 2a #################
library(dplyr)
library(ggbiplot)
library(immunarch)
# PCA for sig SNP
freqs <- read.table("lab.genobaypass_dt_4chr_freq", h=F) #this file is generated with baypass2freqs_cov.py
snpdet <- read.table('lab.snpdet_dt_4chr',h=F)
dat <- cbind(snpdet[,1:2],freqs)
colnames(dat) <- c("BSE-0","BSE-0B","BS3C","BS4C","BSE-1-T1","BSE-1-T2",
				"BSE-1-T2B","BSE-2-T1","BSE-2-T2","BSE-2-T2B","BSE-3-T1","BSE-3-T2","BSE-4-T1","BSE-4-T2","BSE-5-T1",
				"BSE-5-T2","BSE-6-T1","BSE-6-T2","BSE-7-T1","BSE-8-T1","BSE-8-T2","BSE-9-T1","BSE-9-T2","BSE-10-T1",
				"BSE-11-T1","BSE-11-T2","BSE-12-T1","BSE-12-T2")

group <- c(rep("A",2),rep("C",8),"T6","T10","T6","T10","T6","T10","T6","T10","T6","T6","T10","T6","T10","T6","T6","T10","T6","T10")
generation <- as.character(c(rep(0,2),rep(20,2),6,10,10,6,10,10,6,10,6,10,6,10,6,10,6,6,10,6,10,6,6,10,6,10))
sig <- readRDS('lab.sig.RDS')

sig_dat <- merge(dat,sig,by=c("chr","pos"))[,3:30]
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

##number of selected alleles per line
library(ggplot2)
source("parallelism_functions.v2.R")
dat <- readRDS('hap_blocks.rawAFC.RDS')
info <- readRDS('hap_blocks.res.RDS')

res <- calc_num_sel_loci_hb(dat,info,"Neutral999")
gen6 <- unlist(res[[1]])
gen10 <- unlist(res[[2]])
gen10_filled <- c(gen10[1],gen10[2],gen10[3],gen10[4],0,gen10[5],gen10[6],0,gen10[7],gen10[8])

plotdat <- data.frame("pop"=c("Arep3","Brep4","Crep5","Drep6","Erep7","Frep8","Grep9","Hrep10","Irep11","Jrep12"),
                        "generation"=c(rep("Agen6",10),rep("Bgen10",10)),
                        "alleles"=c(gen6,gen10_filled))
p <- ggplot(plotdat, aes(pop,alleles,fill=generation))+
        geom_bar(position="dodge",stat="identity",alpha=0.8) +
        scale_fill_manual(values=c("#7501BC","#C98BF7")) +
        ylim(c(0,150)) +
        geom_hline(yintercept=121, linetype="dashed", color = "red", size=1) +
        theme_classic()
p

### bar plot of Ne per line
library(ggpubr)
plotdat <- data.frame("line"=c("A-BS3C","B-BS4C","C-BSE1","D-BSE2","E-BSE3","F-BSE4","G-BSE5","H-BSE6","I-BSE8","J-BSE9","K-BSE11","L-BSE12"),
                      "Ne"=c(1611.8425,1424.811,2556.4115,2634.097,2271.9895,1873.235,1881.845,1659.2965,1920.6605,1582.5665,1512.775,1502.216),
                    "group"=c(rep("control",4),rep("treatment",8)))

p <- ggdotchart(plotdat, x = "line", y = "Ne",
            color = "group",
            group = "group",
           palette = c("#E7B800","#C98BF7"), # Custom color palette
           sorting = "descending",
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically
           dot.size = 6,                                 # Large dot size
           label = round(plotdat$Ne),                        # Add mpg values as dot labels
           ggtheme = theme_pubr()                        # ggplot2 theme
           )
p
#### Assessing parallelism using pairwise overlap (Jaccard Index)
#
library(ggplot2)
library(dplyr)
library(lme4)
source("extra_scripts/parallelism_functions.R")

# empirical data
# data are formatted to show change from ancestor
dat <- readRDS('full_data_lmm.rawAF.RDS') #this file can be generated from prep_lmm.rawAFC.R
snpdet <- read.table('data/lab.snpdet')
dat$chr <- rep(snpdet[,1],26)
dat$pos <- rep(snpdet[,2],26)

sig <- readRDS('data/lab.sig.RDS')

sig_dat <- merge(dat,sig,by=c("chr","pos"))

####
# jaccard plots under different simulation models
####

info <- readRDS('../haplovalidate_run/hap_blocks.res.RDS')
#empirical data
jaccard_emp <- calc_jaccard_emp(sig_dat,0.01)

#independent SNP sims (sweeps)
# these files are generated using the simulate_poolseq.R script in 'extra_scripts'
jaccard_ind_6 <- c()
jaccard_ind_10 <- c()

ind_sim_files <- list.files(path = ".", pattern = 'freq', full.names=T) #generated from the extra_scripts/simulate_poolseq.R
for (i in 1:length(ind_sim_files)){
	print(i)
	file <- ind_sim_files[i]
	simdat <- read.table(file,h=T)
	res <- calc_jaccard_ind_sims(simdat, info)
	jaccard_ind_6 <- c(jaccard_ind_6,res[[1]])
	jaccard_ind_10 <- c(jaccard_ind_10,res[[2]])
	}

#polygenic shift sims, 30 loci
# these files are generated using the franssen_qt_sims.R script in 'extra_scripts'

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

######
## SNPs
######

# data are formatted to show change from ancestor
# this file can be generated using the prep_lmm.rawAFC.R script
dat <- readRDS('full_data_lmm.rawAF.RDS')
snpdet <- read.table('data/lab.snpdet')
dat$Dovetail_Scaffold <- rep(snpdet[,1],26)
dat$Dovetail_Position <- rep(snpdet[,2],26)

sig <- readRDS('data/lab.sig.RDS')

sig_dat <- merge(dat,sig,by=c("Dovetail_Scaffold","Dovetail_Position"))

## calculate rfs
# creates list where first item is gen 6 and second is gen 10
rfs_emp <- rfs_empirical_snps(sig_dat,0.01)

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
    geom_point(data=sample_data, aes(x=Longitude, y=Latitude, fill=log(MAS)),
               color = "black", shape=21, size=5.0) +
		scale_fill_gradient(low = "#efe1f7", high = "#7501BC") +
    theme(panel.background = element_rect(fill = "#EBF5FB"),
    		axis.title.x=element_blank(),
    		axis.title.y=element_blank())

######Figure 3c #######
### MAF distributions ####
#c. Folded AFS for selected vs non-selected SNP

fold_af <- function(af){
		ifelse(af > 0.5, 1-af,af)
}
res_all <- readRDS('lab.snp_res_dt.RDS')
colnames(res_all)[1:2] <- c("chr","pos")
sig <- readRDS('hap_block_snps.res.RDS')
non_sig <- anti_join(res_all,sig,by=c("chr","pos"))

wild_freqs <- read.table('wild.genobaypass_dt_4chr_freq',h=F)
colnames(wild_freqs) <- c("BB1E","BB2E","GBE","HF1E","HF2E","KIE","RG1E","RG2E","STE")
wild_snpdet <- read.table('wild.snpdet_dt_4chr',h=F)[,1:2]
colnames(wild_snpdet) <- c("chr","pos")
wild_dat <- cbind(wild_snpdet,wild_freqs)

wild_sig <- merge(wild_dat,sig, by=c("chr","pos"))
wild_non_sig <- merge(wild_dat,non_sig, by=c("chr","pos"))

ksres <- ks.test(fold_af(wild_sig$STE), fold_af(wild_non_sig$STE))
ksres$statistic
ksres$p.value

plotdat <- data.frame("Group"=c(rep("Selected",nrow(wild_sig)),rep("Non_Selected",nrow(wild_non_sig))),
					"MAF"= c(fold_af(wild_sig$KIE), fold_af(wild_non_sig$KIE)))


library(plyr); library(dplyr)
library(ggplot2)

cdat <- ddply(plotdat, "Group", summarise, MAF.mean=mean(MAF))
p <- ggplot(plotdat, aes(x = MAF, after_stat(density), color=Group,fill=Group))
p + geom_histogram(binwidth = 0.05, position="dodge")+
	geom_vline(data=cdat, aes(xintercept=MAF.mean,  colour=Group),
               linetype="dashed", size=1) +
    scale_color_manual(values=c("#FFD79D","#C870FF")) +
	scale_fill_manual(values=c("#FFD79D","#C870FF")) +
	theme_classic() +
    theme(legend.position="none", axis.title.x=element_blank(),axis.title.y=element_blank())

ks.test(filter(plotdat, Group=="Selected")$MAF,filter(plotdat, Group=="Non_Selected")$MAF)

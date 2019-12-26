library(lme4)
library(mdscore)
library(glmmTMB)

geno <- read.table("genobaypass_lab_glmm", h=F)

info <- data.frame("Gen"=c(0,6,10,0,6,10,0,6,10,0,6,10,0,6,10,0,6,10,0,6,0,6,10,0,6,10,0,6,0,6,10,0,6,10),
					"Line" = as.factor(c("1","1","1","2","2","2","3","3","3","4","4","4","5","5","5","6","6","6","7","7","8","8","8","9","9","9","10","10","11","11","11","12","12","12")),
					"Treat" = c(rep("Ctrl",6),rep("Sel",28)))				

calc_neff <- function(Nchr, Nrd) {
	neff <- ((Nchr * Nrd) - 1) / (Nchr + Nrd)
	return(neff)
	}

run_lrt_binomial <- function(info, counts){
	dat <- cbind(info, "Acnt"=as.numeric(counts[1,]), "acnt"=as.numeric(counts[2,]))
	rds <- dat$Acnt + dat$acnt
	neffs <- ((100 * rds) - 1) / (100 + rds)
	#neffs <- rep(100,34)
	# only do this on counts larger than 2n
	dat <- cbind(dat, neffs, rds)
	#Acnt_scaled <- c()
	#acnt_scaled <- c()
	#for (x in 1:nrow(dat)) {
	#	if (dat[x,]$rds > dat[x,]$neff){
	#		Acnt_scaled <- c(Acnt_scaled, round((dat[x,]$Acnt / dat[x,]$rds) * dat[x,]$neffs))
	#		acnt_scaled <- c(acnt_scaled, round((dat[x,]$acnt / dat[x,]$rds) * dat[x,]$neffs)) } else {
	#		Acnt_scaled <- c(Acnt_scaled, dat[x,]$Acnt)
	#		acnt_scaled <- c(acnt_scaled, dat[x,]$acnt)
	#		}
	#	}
	#		
	Acnt_scaled <- round((dat$Acnt / dat$rds) * dat$neffs)
	acnt_scaled <- round((dat$acnt / dat$rds) * dat$neffs)
	dat <- cbind(dat, Acnt_scaled, acnt_scaled)
	dat[dat$Acnt_scaled == 0,8] <- 1
	dat[dat$acnt_scaled == 0,9] <- 1
	mod1 <- glmer( cbind(Acnt_scaled, acnt_scaled) ~ 1 + Gen + (Gen+0|Line), family = "binomial", data=dat)
	mod2 <- glmer( cbind(Acnt_scaled, acnt_scaled) ~ Gen*Treat + (Gen+0|Line), family = "binomial", data=dat)
	res <- anova(mod1,mod2)
	return(res$Chisq[2])
	#res <- lr.test(mod1,mod2)
	#return(res$pvalue)
	}

run_lrt_binomial_nobase <- function(info, counts){
	dat <- cbind(info, "Acnt"=as.numeric(counts[1,]), "acnt"=as.numeric(counts[2,]))
	dat <- dat[which(dat$Gen != 0),]
	rds <- dat$Acnt + dat$acnt
	neffs <- ((100 * rds) - 1) / (100 + rds)
	dat <- cbind(dat, neffs, rds)
	Acnt_scaled <- round((dat$Acnt / dat$rds) * dat$neffs)
	acnt_scaled <- round((dat$acnt / dat$rds) * dat$neffs)
	dat <- cbind(dat, Acnt_scaled, acnt_scaled)
	dat[dat$Acnt_scaled == 0,8] <- 1
	dat[dat$acnt_scaled == 0,9] <- 1
	mod1 <- glmer( cbind(Acnt_scaled, acnt_scaled) ~ 1 + Gen + (Gen+0|Line), family = "binomial", data=dat)
	mod2 <- glmer( cbind(Acnt_scaled, acnt_scaled) ~ Gen*Treat + (Gen+0|Line), family = "binomial", data=dat)
	res <- anova(mod1,mod2)
	return(res$Chisq[2])
	}


run_lrt_beta <- function(info, counts){
	dat <- cbind(info, "Acnt"=as.numeric(counts[1,]), "acnt"=as.numeric(counts[2,]))
	#dat <- dat[which(dat$Gen != 0),]
	#dat[dat$Acnt == 0,4] <- 1
	#dat[dat$acnt == 0,5] <- 1
	dat$rds <- dat$Acnt + dat$acnt
	dat$AF <- dat$acnt / dat$rds
	dat[which(dat$AF == 0),7] <- 0.001
	dat[which(dat$AF == 1),7] <- 0.999
	mod1 <- glmmTMB(AF ~ 1 + Gen + (Gen+0|Line), data=dat, beta_family(link="logit"))
	mod2 <- glmmTMB(AF ~ 1+ Gen*Treat + (Gen+0|Line), data=dat, beta_family(link="logit"))
	res <- anova(mod1,mod2)
	return(res$Chisq[2])
	#return(res$'Pr(>F)'[2])
	}

run_lrt_nobase <- function(info, counts){
	dat <- cbind(info, "Acnt"=as.numeric(counts[1,]), "acnt"=as.numeric(counts[2,]))
	dat <- dat[which(dat$Gen != 0),]
	#dat[dat$Acnt == 0,4] <- 1
	#dat[dat$acnt == 0,5] <- 1
	dat$rds <- dat$Acnt + dat$acnt
	dat$AF <- dat$acnt / dat$rds
	dat[which(dat$AF == 0),7] <- 0.001
	dat[which(dat$AF == 1),7] <- 0.999
	mod1 <- glmmTMB(AF ~ 1 + Gen + (Gen+0|Line), data=dat, beta_family(link="logit"))
	mod2 <- glmmTMB(AF ~ 1+ Gen*Treat + (Gen+0|Line), data=dat, beta_family(link="logit"))
	#mod1 <- lm(AF ~ 1 + Gen, data=dat)
	#mod2 <- lm(AF ~ 1+ Gen*Treat, data=dat)
	res <- anova(mod1,mod2)
	return(res$Chisq[2])
	#return(res$'Pr(>F)'[2])
	}

lrt_binomial <- c()
for (i in seq(1,nrow(geno),2)){
	print(i)
	counts <- geno[i:(i+1),]
	res <- run_lrt_binomial(info=info, counts=counts)
	lrt_binomial <- c(lrt_binomial, res)
	}

lrt_binomial_nobase <- c()
for (i in seq(1,nrow(geno),2)){
	print(i)
	counts <- geno[i:(i+1),]
	res <- run_lrt_binomial_nobase(info=info, counts=counts)
	lrt_binomial_nobase <- c(lrt_binomial_nobase, res)
	}


lrt_stat <- c()	
for (i in seq(1,nrow(geno),2)){
	print(i)
	counts <- geno[i:(i+1),]
	res <- run_lrt_beta(info=info, counts=counts)
	lrt_stat <- c(lrt_stat, res)
	}


lrt_nobase_stat <- c()	
for (i in seq(1,nrow(geno),2)){
	print(i)
	counts <- geno[i:(i+1),]
	res <- run_lrt_nobase(info=info, counts=counts)
	lrt_nobase_stat <- c(lrt_nobase_stat, res)
	}


#### Simulations ####
library(poolSeq)

baseAF <- read.table("genobaypass_lab_freq",h=F)[,1]
Cov <- read.table("genobaypass_lab_cov",h=F)
p0=sample(baseAF,10000)

sim_AF_mat <- data.frame(matrix(NA, nrow = 10000, ncol = 36))

sim_Cov_mat <- data.frame(matrix(NA, nrow = 10000, ncol = 36))

for (i in 1:12){
	simTraj <- wf.traj(p0=p0, Ne=900, t=c(0, 6, 10))
	simTraj <- matrix(sample.alleles(simTraj, size=50, mode="individuals", Ncensus=500, ploidy = 2), nrow=nrow(simTraj), dimnames=dimnames(simTraj))
	af <- sample.alleles(simTraj, size=mean(Cov[,1]), mode="coverage")
	simTraj <- matrix(af$p.smpld, nrow=nrow(simTraj), dimnames=dimnames(simTraj))
	sim_AF_mat[,((i*3)-2):(i*3)] <- simTraj
	simCov <- matrix(af$size, nrow=nrow(simTraj), dimnames=dimnames(simTraj))
	sim_Cov_mat[,((i*3)-2):(i*3)] <- simCov
	}
	
## remove BSE7-T2,BSE9-T2
#sim_AF_mat <- sim_AF_mat[,c(1:20,22:29,31:36)]
#sim_cov_mat <- sim_cov_mat[,c(1:20,22:29,31:36)]

## prep for glmm (place BSE0 as ancestor to all others)
sim_AF_mat <- sim_AF_mat[,c(1,2,3,1,5,6,1,8,9,1,11,12,1,14,15,1,17,18,1,20,1,22,23,1,26,27,1,29,1,32,33,1,35,36)]
sim_Cov_mat <- sim_Cov_mat[,c(1,2,3,1,5,6,1,8,9,1,11,12,1,14,15,1,17,18,1,20,1,22,23,1,26,27,1,29,1,32,33,1,35,36)]
write.table(sim_AF_mat,"sim10k_glmm_freq",quote=F,row.names=F)
write.table(sim_Cov_mat,"sim10k_glmm_cov",quote=F,row.names=F)

## prep for LRT12 (only BSE0 and T2)
#sim_AF_mat_lrt <- sim_AF_mat[,c(1,9,1,12,1,15,1,18,1,23,1,27,1,33,1,36)]
#sim_cov_mat_lrt <- sim_Cov_mat[,c(1,9,1,12,1,15,1,18,1,23,1,27,1,33,1,36)]
#write.table(sim_AF_mat_lrt,"sim10k_lrt_freq",quote=F,row.names=F)
#write.table(sim_cov_mat_lrt,"sim10k_lrt_cov",quote=F,row.names=F)

run_lrt_beta <- function(info, counts){
	dat <- cbind(info, "AF"=as.numeric(afs))
	dat[which(dat$AF == 0),4] <- 0.001
	dat[which(dat$AF == 1),4] <- 0.999
	#dat <- dat[which(dat$Gen != 0),]
	mod1 <- glmmTMB(AF ~ 1 + Gen + (Gen+0|Line), data=dat, beta_family(link="logit"))
	mod2 <- glmmTMB(AF ~ 1+ Gen*Treat + (Gen+0|Line), data=dat, beta_family(link="logit"))
	res <- anova(mod1,mod2)
	return(res$Chisq[2])
	#return(res$'Pr(>F)'[2])
	}
	
run_lrt_binomial <- function(info, AF,cov){
	acnt <- round(AF*cov)
	Acnt <- round(cov-(AF*cov))
	neffs <- ((100 * cov) - 1) / (100 + cov)
	dat <- cbind(info, "Acnt"=as.numeric(Acnt), "acnt"=as.numeric(acnt),"rds"=as.numeric(cov),"neffs"=as.numeric(neffs))
	Acnt_scaled <- round((dat$Acnt / dat$rds) * dat$neffs)
	acnt_scaled <- round((dat$acnt / dat$rds) * dat$neffs)
	dat <- cbind(dat, Acnt_scaled, acnt_scaled)
	dat[dat$Acnt_scaled == 0,8] <- 1
	dat[dat$acnt_scaled == 0,9] <- 1
	mod1 <- glmer( cbind(Acnt_scaled, acnt_scaled) ~ 1 + Gen + (Gen+0|Line), family = "binomial", data=dat)
	mod2 <- glmer( cbind(Acnt_scaled, acnt_scaled) ~ Gen*Treat + (Gen+0|Line), family = "binomial", data=dat)
	res <- anova(mod1,mod2)
	return(res$Chisq[2])
	}

run_lrt_binomial_nobase <- function(info, AF,cov){
	acnt <- round(AF*cov)
	Acnt <- round(cov-(AF*cov))
	neffs <- ((100 * cov) - 1) / (100 + cov)
	dat <- cbind(info, "Acnt"=as.numeric(Acnt), "acnt"=as.numeric(acnt),"rds"=as.numeric(cov),"neffs"=as.numeric(neffs))
	dat <- dat[which(dat$Gen != 0),]
	Acnt_scaled <- round((dat$Acnt / dat$rds) * dat$neffs)
	acnt_scaled <- round((dat$acnt / dat$rds) * dat$neffs)
	dat <- cbind(dat, Acnt_scaled, acnt_scaled)
	dat[dat$Acnt_scaled == 0,8] <- 1
	dat[dat$acnt_scaled == 0,9] <- 1
	mod1 <- glmer( cbind(Acnt_scaled, acnt_scaled) ~ 1 + Gen + (Gen+0|Line), family = "binomial", data=dat)
	mod2 <- glmer( cbind(Acnt_scaled, acnt_scaled) ~ Gen*Treat + (Gen+0|Line), family = "binomial", data=dat)
	res <- anova(mod1,mod2)
	return(res$Chisq[2])
	}


run_lrt_nobase_beta <- function(info, counts){
	dat <- cbind(info, "AF"=as.numeric(afs))
	dat[which(dat$AF == 0),4] <- 0.001
	dat[which(dat$AF == 1),4] <- 0.999
	dat <- dat[which(dat$Gen != 0),]
	mod1 <- glmmTMB(AF ~ 1 + Gen + (Gen+0|Line), data=dat, beta_family(link="logit"))
	mod2 <- glmmTMB(AF ~ 1+ Gen*Treat + (Gen+0|Line), data=dat, beta_family(link="logit"))
	res <- anova(mod1,mod2)
	return(res$Chisq[2])
	#return(res$'Pr(>F)'[2])
	}


lrt_stat <- c()	
for (i in 1:nrow(sim_AF_mat)){
	print(i)
	afs <- sim_AF_mat[i,]
	res <- run_lrt_beta(info=info, counts=afs)
	lrt_stat <- c(lrt_stat, res)
	}

lrt_nobase_stat <- c()	
for (i in 1:nrow(sim_AF_mat)){
	print(i)
	afs <- sim_AF_mat[i,]
	res <- run_lrt_nobase_beta(info=info, counts=afs)
	lrt_nobase_stat <- c(lrt_nobase_stat, res)
	}

binomial_lrt_stat <- c()	
for (i in 1:nrow(sim_AF_mat)){
	print(i)
	afs <- sim_AF_mat[i,]
	cov <- sim_Cov_mat[i,] 
	res <- run_lrt_binomial(info=info, AF=afs,cov=cov)
	binomial_lrt_stat <- c(binomial_lrt_stat, res)
	}

binomial_lrt_nobase <- c()	
for (i in 1:nrow(sim_AF_mat)){
	print(i)
	afs <- sim_AF_mat[i,]
	cov <- sim_Cov_mat[i,] 
	res <- run_lrt_binomial_nobase(info=info, AF=afs,cov=cov)
	binomial_lrt_nobase <- c(binomial_lrt_nobase, res)
	}

#########################################################################################
###### 
###### Using new controls: BSE0B, BSE3C, BSE4C, BSE-1-T2B, BSE-2-T2B
######
#########################################################################################


## Simple binomial GLM where all 'non-selected' / control lines are 'Gen0' and the rest are gen 6 and 10
#  Take into account random 'line' effect just for gen 6 and 10 beakers with two samples 


library(lme4)

geno <- read.table("genobaypass_lab_glmm", h=F)

info <- data.frame("Gen"=c(0,0,0,0,0,0,0,0,0,0,6,10,6,10,6,10,6,10,6,6,10,6,10,6,6,10,6,10),
					"Batch" = as.factor(c("1","3","3","3","2","1","3","2","1","3","2","1","2","1","2","1","2","1","2","2","1","2","1","2","2","1","2","1")),
					"Line" = as.factor(c(rep(NA,10),"3","3","4","4","5","5","6","6",NA,"8","8","9","9",NA,"11","11","12","12")))

run_lrt_binomial <- function(info, counts){
	dat <- cbind(info, "Acnt"=as.numeric(counts[1,]), "acnt"=as.numeric(counts[2,]))
	rds <- dat$Acnt + dat$acnt
	neffs <- ((100 * rds) - 1) / (100 + rds)
	dat <- cbind(dat, neffs, rds)
	Acnt_scaled <- round((dat$Acnt / dat$rds) * dat$neffs)
	acnt_scaled <- round((dat$acnt / dat$rds) * dat$neffs)
	dat <- cbind(dat, Acnt_scaled, acnt_scaled)
	dat[dat$Acnt_scaled == 0,8] <- 1
	dat[dat$acnt_scaled == 0,9] <- 1
	## Assume random slope and / or intercept for each line
	mod <- glm( cbind(Acnt_scaled, acnt_scaled) ~ Gen, family = "binomial", data=dat)
	res <- coef(summary(mod))[2,4]
	return(res)
	}

binomial_stat <- c()
for (i in seq(1,nrow(geno),2)){
	print(i)
	counts <- geno[i:(i+1),]
	res <- run_lrt_binomial(info=info, counts=counts)
	binomial_stat <- c(binomial_stat, res)
	}

###################################
### Use BSE0 and BSE0B to simulate base for everyone
###################################

library(lme4)
library(lmtest)

geno <- read.table("genobaypass_lab_glmm", h=F)

info <- data.frame("Gen"=c(0,20,0,20,0,6,10,10,0,6,10,10,0,6,10,0,6,10,0,6,10,0,6,10,0,6,0,6,10,0,6,10,0,6,0,6,10,0,6,10),
					"NChrom"=c(200,100,200,100,200,100,100,100,200,100,100,100,200,100,100,200,100,100,200,100,100,200,100,100,200,100,200,100,100,200,100,100,200,100,200,100,100,200,100,100),
					"Line" = as.factor(c("BS3C","BS3C","BS4C","BS4C","1","1","1","1","2","2","2","2","3","3","3","4","4","4","5","5","5","6","6","6","7","7","8","8","8","9","9","9","10","10","11","11","11","12","12","12")),
					"Treat" = c(rep("Ctrl",12),rep("Sel",28)))	
					
								
run_lrt_binomial <- function(counts,info){
	base_Acnt <- counts[1,1]+counts[1,2]
	base_acnt <- counts[2,1]+counts[2,2]
	## simulate 'base' data for the 14 samples based on the base data
	sim_Acnt <- rbinom(14,round((base_Acnt + base_acnt)/2), (base_Acnt / (base_Acnt + base_acnt)))
	sim_acnt <- rbinom(14,round((base_Acnt + base_acnt)/2), (base_acnt / (base_Acnt + base_acnt)))
	dat <- cbind(info, 
			"Acnt"=as.numeric(c(sim_Acnt[1],counts[1,3],sim_Acnt[2],counts[1,4],sim_Acnt[3],counts[1,5:7],
			sim_Acnt[4],counts[1,8:10],sim_Acnt[5],counts[1,11:12],sim_Acnt[6],counts[1,13:14],
			sim_Acnt[7],counts[1,15:16],sim_Acnt[8],counts[1,17:18],sim_Acnt[9],counts[1,19],
			sim_Acnt[10],counts[1,20:21],sim_Acnt[11],counts[1,22:23],sim_Acnt[12],counts[1,24],
			sim_Acnt[13],counts[1,25:26],sim_Acnt[14],counts[1,27:28])), 
			"acnt"=as.numeric(c(sim_acnt[1],counts[2,3],sim_acnt[2],counts[2,4],sim_acnt[3],counts[2,5:7],
			sim_acnt[4],counts[2,8:10],sim_acnt[5],counts[2,11:12],sim_acnt[6],counts[2,13:14],
			sim_acnt[7],counts[2,15:16],sim_acnt[8],counts[2,17:18],sim_acnt[9],counts[2,19],
			sim_acnt[10],counts[2,20:21],sim_acnt[11],counts[2,22:23],sim_acnt[12],counts[2,24],
			sim_acnt[13],counts[2,25:26],sim_acnt[14],counts[2,27:28])))
	rds <- dat$Acnt + dat$acnt
	#neffs <- ((dat$NChrom * rds) - 1) / (dat$NChrom + rds)
	neffs <- ((100 * rds) - 1) / (100 + rds)
	dat <- cbind(dat, neffs, rds)
	Acnt_scaled <- round((dat$Acnt / dat$rds) * dat$neffs)
	acnt_scaled <- round((dat$acnt / dat$rds) * dat$neffs)
	dat <- cbind(dat, Acnt_scaled, acnt_scaled)
	dat[dat$Acnt_scaled == 0,9] <- 1
	dat[dat$acnt_scaled == 0,10] <- 1
	## assume random slope, same intercept for each line
	#mod1 <- glmer( cbind(Acnt_scaled, acnt_scaled) ~ Gen + (Gen+0|Line), family = "binomial", data=dat)
	#mod2 <- glmer( cbind(Acnt_scaled, acnt_scaled) ~ Gen*Treat + (Gen+0|Line), family = "binomial", data=dat)
	#res <- anova(mod1,mod2)
	#return(res$Chisq[2])
	#### Dont take random line effect into account
	mod1 <- glm( cbind(Acnt_scaled, acnt_scaled) ~ Gen, family = "binomial", data=dat)
	mod2 <- glm( cbind(Acnt_scaled, acnt_scaled) ~ Gen*Treat, family = "binomial", data=dat)
	res <- lrtest(mod1,mod2)
	return(res$Chisq[2])
	}

lrt_binomial <- c()
for (i in seq(1,nrow(geno),2)){
	print(i)
	counts <- geno[i:(i+1),]
	res <- run_lrt_binomial(info=info, counts=counts)
	lrt_binomial <- c(lrt_binomial, res)
	}
	
pvals <- sapply(lrt_binomial, pchisq, df=2)
pvals <- 1-pvals

############### Simulations ##########
library(poolSeq)

baseAF <- read.table("genobaypass_lab_freq",h=F)[,1]
Cov <- read.table("genobaypass_lab_cov",h=F)
p0=sample(baseAF,10000)

sim_AF_mat <- data.frame(matrix(NA, nrow = 10000, ncol = 56))

sim_Cov_mat <- data.frame(matrix(NA, nrow = 10000, ncol = 56))

for (i in 1:14){
	simTraj <- wf.traj(p0=p0, Ne=900, t=c(0, 6, 10, 20))
	simTraj <- matrix(sample.alleles(simTraj, size=50, mode="individuals", Ncensus=500, ploidy = 2), nrow=nrow(simTraj), dimnames=dimnames(simTraj))
	af <- sample.alleles(simTraj, size=mean(Cov[,1]), mode="coverage")
	simTraj <- matrix(af$p.smpld, nrow=nrow(simTraj), dimnames=dimnames(simTraj))
	sim_AF_mat[,((i*4)-3):(i*4)] <- simTraj
	simCov <- matrix(af$size, nrow=nrow(simTraj), dimnames=dimnames(simTraj))
	sim_Cov_mat[,((i*4)-3):(i*4)] <- simCov
	}
											
sim_AF_mat_glmm_imputeBase <- sim_AF_mat[,c(1,4,5,8,9,10,11,11,13,14,15,15,17,18,19,21,22,23,25,26,27,29,30,31,33,34,37,38,39,41,42,43,45,46,49,50,51,53,54,55)]
sim_Cov_mat_glmm_imputeBase <- sim_Cov_mat[,c(1,4,5,8,9,10,11,11,13,14,15,15,17,18,19,21,22,23,25,26,27,29,30,31,33,34,37,38,39,41,42,43,45,46,49,50,51,53,54,55)]
write.table(sim_AF_mat,"sim10k_glmm_imputeBase_freq",quote=F,row.names=F)
write.table(sim_Cov_mat,"sim10k_glmm_imputeBase_cov",quote=F,row.names=F)

sim_AF_mat_glm_gen <- sim_AF_mat[,c(1,5,4,8,10,11,11,14,15,15,18,19,22,23,26,27,30,31,34,38,39,42,43,46,50,51,54,55)]
sim_Cov_mat_glm_gen <- sim_Cov_mat[,c(1,5,4,8,10,11,11,14,15,15,18,19,22,23,26,27,30,31,34,38,39,42,43,46,50,51,54,55)]
write.table(sim_AF_mat,"sim10k_glm_gen_freq",quote=F,row.names=F)
write.table(sim_Cov_mat,"sim10k_glm_gen_cov",quote=F,row.names=F)

### binomial GLM LRT with imputed base
info <- data.frame("Gen"=c(0,20,0,20,0,6,10,10,0,6,10,10,0,6,10,0,6,10,0,6,10,0,6,10,0,6,0,6,10,0,6,10,0,6,0,6,10,0,6,10),
					"Line" = as.factor(c("BS3C","BS3C","BS4C","BS4C","1","1","1","1","2","2","2","2","3","3","3","4","4","4","5","5","5","6","6","6","7","7","8","8","8","9","9","9","10","10","11","11","11","12","12","12")),
					"Treat" = c(rep("Ctrl",10),rep("Sel",18)))	

run_lrt_binomial <- function(info, AF,cov){
	acnt <- round(AF*cov)
	Acnt <- round(cov-(AF*cov))
	neffs <- ((100 * cov) - 1) / (100 + cov)
	dat <- cbind(info, "Acnt"=as.numeric(Acnt), "acnt"=as.numeric(acnt),"rds"=as.numeric(cov),"neffs"=as.numeric(neffs))
	Acnt_scaled <- round((dat$Acnt / dat$rds) * dat$neffs)
	acnt_scaled <- round((dat$acnt / dat$rds) * dat$neffs)
	dat <- cbind(dat, Acnt_scaled, acnt_scaled)
	dat[dat$Acnt_scaled == 0,8] <- 1
	dat[dat$acnt_scaled == 0,9] <- 1
	mod1 <- glmer( cbind(Acnt_scaled, acnt_scaled) ~ 1 + Gen + (Gen+0|Line), family = "binomial", data=dat)
	mod2 <- glmer( cbind(Acnt_scaled, acnt_scaled) ~ Gen*Treat + (Gen+0|Line), family = "binomial", data=dat)
	res <- anova(mod1,mod2)
	return(res$Chisq[2])
	}

sim_AF_mat_glmm_imputeBase <- sim_AF_mat_glmm_imputeBase[which(rowSums(sim_AF_mat_glmm_imputeBase) > 0),]
sim_Cov_mat_glmm_imputeBase <- sim_Cov_mat_glmm_imputeBase[which(rowSums(sim_AF_mat_glmm_imputeBase) > 0),]


binomial_lrt_stat <- c()	
for (i in 1:nrow(sim_AF_mat_glmm_imputeBase)){
	print(i)
	afs <- sim_AF_mat_glmm_imputeBase[i,]
	cov <- sim_Cov_mat_glmm_imputeBase[i,] 
	res <- try(run_lrt_binomial(info=info, AF=afs,cov=cov), silent=T)
	binomial_lrt_stat <- c(binomial_lrt_stat, res)
	}
lrt_stat <- sapply(binomial_lrt_stat, as.numeric)
pvals <- sapply(lrt_stat, pchisq, df=2)
pvals <- 1-pvals

### binomial GLM LRT without imputed base
library(lme4)
sim_AF <- read.table('sim10k_glm_gen_freq',h=F)
sim_Cov <- read.table('sim10k_glm_gen_cov',h=F)

info <- data.frame("Gen"=c(0,20,0,20,0,6,10,10,0,6,10,10,0,6,10,0,6,10,0,6,10,0,6,10,0,6,0,6,10,0,6,10,0,6,0,6,10,0,6,10),
					"Line" = as.factor(c("BS3C","BS3C","BS4C","BS4C","1","1","1","1","2","2","2","2","3","3","3","4","4","4","5","5","5","6","6","6","7","7","8","8","8","9","9","9","10","10","11","11","11","12","12","12")),
					"Treat" = c(rep("Ctrl",12),rep("Sel",28)))	

run_lrt_binomial <- function(info, AF,cov){
	acnt <- round(AF*cov)
	Acnt <- round(cov-(AF*cov))
	neffs <- ((100 * cov) - 1) / (100 + cov)
	dat <- cbind(info, "Acnt"=as.numeric(Acnt), "acnt"=as.numeric(acnt),"rds"=as.numeric(cov),"neffs"=as.numeric(neffs))
	Acnt_scaled <- round((dat$Acnt / dat$rds) * dat$neffs)
	acnt_scaled <- round((dat$acnt / dat$rds) * dat$neffs)
	dat <- cbind(dat, Acnt_scaled, acnt_scaled)
	dat[dat$Acnt_scaled == 0,8] <- 1
	dat[dat$acnt_scaled == 0,9] <- 1
	mod1 <- glmer( cbind(Acnt_scaled, acnt_scaled) ~ 1 + Gen + (Gen+0|Line), family = "binomial", data=dat)
	mod2 <- glmer( cbind(Acnt_scaled, acnt_scaled) ~ Gen*Treat + (Gen+0|Line), family = "binomial", data=dat)
	res <- anova(mod1,mod2)
	return(res$Chisq[2])
	}

sim_AF_mat_glmm_imputeBase <- sim_AF_mat_glmm_imputeBase[which(rowSums(sim_AF_mat_glmm_imputeBase) > 0),]
sim_Cov_mat_glmm_imputeBase <- sim_Cov_mat_glmm_imputeBase[which(rowSums(sim_AF_mat_glmm_imputeBase) > 0),]


binomial_lrt_stat <- c()	
for (i in 1:nrow(sim_AF_mat_glmm_imputeBase)){
	print(i)
	afs <- sim_AF_mat_glmm_imputeBase[i,]
	cov <- sim_Cov_mat_glmm_imputeBase[i,] 
	res <- try(run_lrt_binomial(info=info, AF=afs,cov=cov), silent=T)
	binomial_lrt_stat <- c(binomial_lrt_stat, res)
	}
lrt_stat <- sapply(binomial_lrt_stat, as.numeric)
pvals <- sapply(lrt_stat, pchisq, df=2)
pvals <- 1-pvals

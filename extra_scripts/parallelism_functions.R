# functions to calculate Jaccard statistics and the replicate frequency spectrum (sensu Barghi et al. 2019)
# for simulated and real data from Stern, Diaz, and Lee 2020
# Mostly for organizational purposes. Not great for general use

library(lme4)
library(dplyr)

calc_jaccard_ind_sims <- function(sig_simdat,afc){
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
	  	minafc <- afc
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

		gen6_names <- c("BSE3_6","BSE4_6","BSE5_6","BSE6_6","BSE7_6","BSE8_6","BSE9_6","BSE10_6","BSE11_6","BSE12_6")
		gen10_names <- c("BSE3_10","BSE4_10","BSE5_10","BSE6_10","BSE7_10","BSE8_10","BSE9_10","BSE10_10","BSE11_10","BSE12_10")
		gen6_pops <- list(BSE3_6,BSE4_6,BSE5_6,BSE6_6,BSE7_6,BSE8_6,BSE9_6,BSE10_6,BSE11_6,BSE12_6)
		names(gen6_pops) <- gen6_names
		gen10_pops <- list(BSE3_10,BSE4_10,BSE5_10,BSE6_10,BSE7_10,BSE8_10,BSE9_10,BSE10_10,BSE11_10,BSE12_10)
		names(gen10_pops) <- gen10_names

		gen6_jac <- matrix(nrow = length(gen6_pops), ncol = length(gen6_pops),
								dimnames = list(c(gen6_names),c(gen6_names)))

		gen10_jac <- matrix(nrow = length(gen10_pops), ncol = length(gen10_pops),
								dimnames = list(c(gen10_names),c(gen10_names)))

		for (i in 1:10){
			for (x in 1:10){
				if (i != x){
					pop1 <- gen6_names[[i]]
					pop2 <- gen6_names[[x]]
					jac <- length(intersect(gen6_pops[[pop1]],gen6_pops[[pop2]])) / length(union(gen6_pops[[pop1]],gen6_pops[[pop2]]))
					#gen6_jac <- c(gen6_jac,jac)
					gen6_jac[pop1,pop2] <- jac
					}
				}
			}


		for (i in 1:10){
				for (x in 1:10){
					if (i != x){
						pop1 <- gen10_names[[i]]
						pop2 <- gen10_names[[x]]
						jac <- length(intersect(gen10_pops[[pop1]],gen10_pops[[pop2]])) / length(union(gen10_pops[[pop1]],gen10_pops[[pop2]]))
						#gen10_jac <- c(gen10_jac,jac)
						gen10_jac[pop1,pop2] <- jac
						}
					}
				}

		return(list(gen6_jac[lower.tri(gen6_jac, diag = FALSE)],gen10_jac[lower.tri(gen10_jac, diag = FALSE)]))
}

calc_jaccard_qt_sims <- function(simdat,afc){
	Rep1_6 <- c()
	Rep1_10 <- c()
	Rep2_6 <- c()
	Rep2_10 <- c()
	Rep3_6 <- c()
	Rep3_10 <- c()
	Rep4_6 <- c()
	Rep4_10 <- c()
	Rep5_6 <- c()
	Rep5_10 <- c()
	Rep6_6 <- c()
	Rep6_10 <- c()
	Rep7_6 <- c()
	Rep7_10 <- c()
	Rep8_6 <- c()
	Rep8_10 <- c()
	Rep9_6 <- c()
	Rep9_10 <- c()
	Rep10_6 <- c()
	Rep10_10 <- c()


	for (i in 1:max(simdat$snp_id)){
		#print(i)
		dt <- filter(simdat,snp_id==i)
		anc <- filter(dt, replicate==1 & generation==0)$frequency
		gen6 <- filter(dt, generation==6)$frequency - anc
		gen10 <- filter(dt, generation==10)$frequency - anc
	  	minafc <- afc
		##gen6
		if (gen6[1]>minafc){
			Rep1_6 <- c(Rep1_6,dt$snp_id[1])
			}
		if (gen6[2]>minafc){
			Rep2_6 <- c(Rep2_6,dt$snp_id[1])
			}
		if (gen6[3]>minafc){
			Rep3_6 <- c(Rep3_6,dt$snp_id[1])
			}
		if (gen6[4]>minafc){
			Rep4_6 <- c(Rep4_6,dt$snp_id[1])
			}
		if (gen6[5]>minafc){
			Rep5_6 <- c(Rep5_6,dt$snp_id[1])
			}
		if (gen6[6]>minafc){
			Rep6_6 <- c(Rep6_6,dt$snp_id[1])
			}
		if (gen6[7]>minafc){
			Rep7_6 <- c(Rep7_6,dt$snp_id[1])
			}
		if (gen6[8]>minafc){
			Rep8_6 <- c(Rep8_6,dt$snp_id[1])
			}
		if (gen6[9]>minafc){
			Rep9_6 <- c(Rep9_6,dt$snp_id[1])
			}
		if (gen6[10]>minafc){
			Rep10_6 <- c(Rep10_6,dt$snp_id[1])
			}
		##gen10
		if (gen10[1] >minafc){
			Rep1_10 <- c(Rep1_10,dt$snp_id[1])
			}
		if (gen10[2] >minafc){
			Rep2_10 <- c(Rep2_10,dt$snp_id[1])
			}
		if (gen10[3] >minafc){
			Rep3_10 <- c(Rep3_10,dt$snp_id[1])
			}
		if (gen10[4] >minafc){
			Rep4_10 <- c(Rep4_10,dt$snp_id[1])
			}
		if (gen10[5] >minafc){
			Rep5_10 <- c(Rep5_10,dt$snp_id[1])
			}
		if (gen10[6] >minafc){
			Rep6_10 <- c(Rep6_10,dt$snp_id[1])
			}
		if (gen10[7] >minafc){
			Rep7_10 <- c(Rep7_10,dt$snp_id[1])
			}
		if (gen10[8] >minafc){
			Rep8_10 <- c(Rep8_10,dt$snp_id[1])
			}
		if (gen10[9] >minafc){
			Rep9_10 <- c(Rep9_10,dt$snp_id[1])
			}
		if (gen10[10] >minafc){
			Rep10_10 <- c(Rep10_10,dt$snp_id[1])
			}
		}

		gen6_names <- c("Rep1_6","Rep2_6","Rep3_6","Rep4_6","Rep5_6","Rep6_6","Rep7_6","Rep8_6","Rep9_6","Rep10_6")
		gen10_names <- c("Rep1_10","Rep2_10","Rep3_10","Rep4_10","Rep5_10","Rep6_10","Rep7_10","Rep8_10","Rep9_10","Rep10_10")
		gen6_pops <- list(Rep1_6,Rep2_6,Rep3_6,Rep4_6,Rep5_6,Rep6_6,Rep7_6,Rep8_6,Rep9_6,Rep10_6)
		names(gen6_pops) <- gen6_names
		gen10_pops <- list(Rep1_10,Rep2_10,Rep3_10,Rep4_10,Rep5_10,Rep6_10,Rep7_10,Rep8_10,Rep9_10,Rep10_10)
		names(gen10_pops) <- gen10_names

		gen6_jac <- matrix(nrow = length(gen6_pops), ncol = length(gen6_pops),
								dimnames = list(c(gen6_names),c(gen6_names)))

		gen10_jac <- matrix(nrow = length(gen10_pops), ncol = length(gen10_pops),
								dimnames = list(c(gen10_names),c(gen10_names)))

		for (i in 1:10){
			for (x in 1:10){
				if (i != x){
					pop1 <- gen6_names[[i]]
					pop2 <- gen6_names[[x]]
					jac <- length(intersect(gen6_pops[[pop1]],gen6_pops[[pop2]])) / length(union(gen6_pops[[pop1]],gen6_pops[[pop2]]))
					#gen6_jac <- c(gen6_jac,jac)
					gen6_jac[pop1,pop2] <- jac
					}
				}
			}


		for (i in 1:10){
				for (x in 1:10){
					if (i != x){
						pop1 <- gen10_names[[i]]
						pop2 <- gen10_names[[x]]
						jac <- length(intersect(gen10_pops[[pop1]],gen10_pops[[pop2]])) / length(union(gen10_pops[[pop1]],gen10_pops[[pop2]]))
						#gen10_jac <- c(gen10_jac,jac)
						gen10_jac[pop1,pop2] <- jac
						}
					}
				}

		return(list(gen6_jac[lower.tri(gen6_jac, diag = FALSE)],gen10_jac[lower.tri(gen10_jac, diag = FALSE)]))

}

calc_jaccard_emp <- function(sig_dat,afc){
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
		meanafc <- dt$selCoef[[1]]
		minafc <- afc
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

	gen6_names <- c("BSE3_6","BSE4_6","BSE5_6","BSE6_6","BSE7_6","BSE8_6","BSE9_6","BSE10_6","BSE11_6","BSE12_6")
	gen10_names <- c("BSE3_10","BSE4_10","BSE5_10","BSE6_10","BSE8_10","BSE9_10","BSE11_10","BSE12_10")
	gen6_pops <- list(BSE3_6,BSE4_6,BSE5_6,BSE6_6,BSE7_6,BSE8_6,BSE9_6,BSE10_6,BSE11_6,BSE12_6)
	names(gen6_pops) <- gen6_names
	gen10_pops <- list(BSE3_10,BSE4_10,BSE5_10,BSE6_10,BSE8_10,BSE9_10,BSE11_10,BSE12_10)
	names(gen10_pops) <- gen10_names

	gen6_jac <- matrix(nrow = length(gen6_pops), ncol = length(gen6_pops),
							dimnames = list(c(gen6_names),c(gen6_names)))

	gen10_jac <- matrix(nrow = length(gen10_pops), ncol = length(gen10_pops),
							dimnames = list(c(gen10_names),c(gen10_names)))

	for (i in 1:10){
		for (x in 1:10){
			if (i != x){
				pop1 <- gen6_names[[i]]
				pop2 <- gen6_names[[x]]
				jac <- length(intersect(gen6_pops[[pop1]],gen6_pops[[pop2]])) / length(union(gen6_pops[[pop1]],gen6_pops[[pop2]]))
				#gen6_jac <- c(gen6_jac,jac)
				gen6_jac[pop1,pop2] <- jac
				}
			}
		}


	for (i in 1:8){
			for (x in 1:8){
				if (i != x){
					pop1 <- gen10_names[[i]]
					pop2 <- gen10_names[[x]]
					jac <- length(intersect(gen10_pops[[pop1]],gen10_pops[[pop2]])) / length(union(gen10_pops[[pop1]],gen10_pops[[pop2]]))
					#gen10_jac <- c(gen10_jac,jac)
					gen10_jac[pop1,pop2] <- jac
					}
				}
			}

	return(list(gen6_jac[lower.tri(gen6_jac, diag = FALSE)],gen10_jac[lower.tri(gen10_jac, diag = FALSE)]))
}

rfs_empirical <- function(sig_dat, afc){
	prop_pop6 <- c()
	prop_pop10 <- c()

	for (i in 1:length(unique(sig_dat$SNP))){
		snpnum <- unique(sig_dat$SNP)[i]
		#print(snpnum)
		dt <- filter(sig_dat,SNP==snpnum)
		gen6 <- filter(dt, Treat=="Treatment" & Generation==6)$value
		gen10 <- filter(dt, Treat=="Treatment" & Generation==10)$value
		##gen6
		afchange <- afc
		gen6_pos <- gen6[gen6 >=afchange]
		gen6_neg <- gen6[gen6 <=afchange]
		if (dt$selCoef[[1]] > 0){
			prop_pop6 <- c(prop_pop6,length(gen6_pos)/10)
			} else{
			prop_pop6 <- c(prop_pop6,length(gen6_neg)/10)
			}
		##gen10
		gen10_pos <- gen10[gen10 >=afchange]
		gen10_neg <- gen10[gen10 <=afchange]
		if (dt$selCoef[[1]] > 0){
			prop_pop10 <- c(prop_pop10,length(gen10_pos)/8)
			} else{
			prop_pop10 <- c(prop_pop10,length(gen10_neg)/8)
			}
		}
	return(list(prop_pop6,prop_pop10))
}

rfs_ind_sims <- function(sig_simdat, afc){
	prop_pop6_sim <- c()
	prop_pop10_sim <- c()

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
	  	estS <- coef(sum)[2,1]
	  	afchange <- afc
	  	#gen6
	  	gen6_pos <- gen6[gen6 >=afchange]
		gen6_neg <- gen6[gen6 <=afchange]
		if (estS > 0){
			prop_pop6_sim <- c(prop_pop6_sim,length(gen6_pos)/10)
			} else{
			prop_pop6_sim <- c(prop_pop6_sim,length(gen6_neg)/10)
			}
		##gen10
		gen10_pos <- gen10[gen10 >=afchange]
		gen10_neg <- gen10[gen10 <=afchange]
		if (estS > 0){
			prop_pop10_sim <- c(prop_pop10_sim,length(gen10_pos)/10)
			} else{
			prop_pop10_sim <- c(prop_pop10_sim,length(gen10_neg)/10)
			}
		}
	return(list(prop_pop6_sim,prop_pop10_sim))
}

rfs_qt_sims <- function(sig_simdat, afc){
	prop_pop6_sim <- c()
	prop_pop10_sim <- c()

	for (i in 1:max(sig_simdat$snp_id)){
		#print(i)
		dt <- filter(simdat,snp_id==i)
		anc <- filter(dt, replicate==1 & generation==0)$frequency
		gen6 <- filter(dt, generation==6)$frequency - anc
		gen10 <- filter(dt, generation==10)$frequency - anc
	  afchange <- afc
	 	#gen6
	  gen6_pos <- gen6[gen6 >=afchange]
		prop_pop6_sim <- c(prop_pop6_sim,length(gen6_pos)/10)
		##gen10
		gen10_pos <- gen10[gen10 >=afchange]
		prop_pop10_sim <- c(prop_pop10_sim,length(gen10_pos)/10)
		}

	return(list(prop_pop6_sim,prop_pop10_sim))
}

calc_num_sel_loci <- function(sig_dat,afc){
	#for each replicate, calculate number of loci with frequency change above min afc
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
		meanafc <- dt$selCoef[[1]]
		minafc <- afc
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

	gen6_nloci <- lapply(gen6_pops, length)
	gen10_nloci <- lapply(gen10_pops, length)
	return(list(gen6_nloci,gen10_nloci))
}

jaccard_redundant <- function(nloci){
	#nloci is the list of num loci generated with calc_num_sel_loci
	loc_nums <- 1:1156
	BSE3_6 <- sample(loc_nums,nloci[[1]][[1]])
	BSE4_6 <- sample(loc_nums,nloci[[1]][[2]])
	BSE5_6 <- sample(loc_nums,nloci[[1]][[3]])
	BSE6_6 <- sample(loc_nums,nloci[[1]][[4]])
	BSE7_6 <- sample(loc_nums,nloci[[1]][[5]])
	BSE8_6 <- sample(loc_nums,nloci[[1]][[6]])
	BSE9_6 <- sample(loc_nums,nloci[[1]][[7]])
	BSE10_6 <- sample(loc_nums,nloci[[1]][[8]])
	BSE11_6 <- sample(loc_nums,nloci[[1]][[9]])
	BSE12_6 <- sample(loc_nums,nloci[[1]][[10]])

	BSE3_10 <- sample(loc_nums,nloci[[2]][[1]])
	BSE4_10 <- sample(loc_nums,nloci[[2]][[2]])
	BSE5_10 <- sample(loc_nums,nloci[[2]][[3]])
	BSE6_10 <- sample(loc_nums,nloci[[2]][[4]])
	BSE8_10 <- sample(loc_nums,nloci[[2]][[5]])
	BSE9_10 <- sample(loc_nums,nloci[[2]][[6]])
	BSE11_10 <- sample(loc_nums,nloci[[2]][[7]])
	BSE12_10 <- sample(loc_nums,nloci[[2]][[8]])

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

	return(list(gen6_jac,gen10_jac))

}

rfs_redundant <- function(nloci){
	prop_pop6 <- c()
	prop_pop10 <- c()

	loc_nums <- 1:1156
	BSE3_6 <- sample(loc_nums,nloci[[1]][[1]])
	BSE4_6 <- sample(loc_nums,nloci[[1]][[2]])
	BSE5_6 <- sample(loc_nums,nloci[[1]][[3]])
	BSE6_6 <- sample(loc_nums,nloci[[1]][[4]])
	BSE7_6 <- sample(loc_nums,nloci[[1]][[5]])
	BSE8_6 <- sample(loc_nums,nloci[[1]][[6]])
	BSE9_6 <- sample(loc_nums,nloci[[1]][[7]])
	BSE10_6 <- sample(loc_nums,nloci[[1]][[8]])
	BSE11_6 <- sample(loc_nums,nloci[[1]][[9]])
	BSE12_6 <- sample(loc_nums,nloci[[1]][[10]])

	BSE3_10 <- sample(loc_nums,nloci[[2]][[1]])
	BSE4_10 <- sample(loc_nums,nloci[[2]][[2]])
	BSE5_10 <- sample(loc_nums,nloci[[2]][[3]])
	BSE6_10 <- sample(loc_nums,nloci[[2]][[4]])
	BSE8_10 <- sample(loc_nums,nloci[[2]][[5]])
	BSE9_10 <- sample(loc_nums,nloci[[2]][[6]])
	BSE11_10 <- sample(loc_nums,nloci[[2]][[7]])
	BSE12_10 <- sample(loc_nums,nloci[[2]][[8]])

	gen6_pops <- list(BSE3_6,BSE4_6,BSE5_6,BSE6_6,BSE7_6,BSE8_6,BSE9_6,BSE10_6,BSE11_6,BSE12_6)
	gen10_pops <- list(BSE3_10,BSE4_10,BSE5_10,BSE6_10,BSE8_10,BSE9_10,BSE11_10,BSE12_10)

	for (i in 1:1156){
		gen6 <- unlist(lapply(gen6_pops, function(x) i %in% x))
		prop_pop6 <- c(prop_pop6,length(gen6[gen6])/10)
		gen10 <- unlist(lapply(gen10_pops, function(x) i %in% x))
		prop_pop10 <- c(prop_pop10,length(gen6[gen6])/8)
	}

	return(list(prop_pop6,prop_pop10))
}

af_sync <- function(dat){
	counts <- strsplit(dat,":")
	a <- as.numeric(counts[[1]][1])
	t <- as.numeric(counts[[1]][2])
	af <- a / (a+t)
	return(af)
}

calc_jaccard_mimicree <- function(simdat,afc) {
	Rep1_6 <- c()
	Rep1_10 <- c()
	Rep2_6 <- c()
	Rep2_10 <- c()
	Rep3_6 <- c()
	Rep3_10 <- c()
	Rep4_6 <- c()
	Rep4_10 <- c()
	Rep5_6 <- c()
	Rep5_10 <- c()
	Rep6_6 <- c()
	Rep6_10 <- c()
	Rep7_6 <- c()
	Rep7_10 <- c()
	Rep8_6 <- c()
	Rep8_10 <- c()
	Rep9_6 <- c()
	Rep9_10 <- c()
	Rep10_6 <- c()
	Rep10_10 <- c()

	for (i in 1:nrow(simdat)){
		dat <- filter(simdat,V1==i)
		dat <- t(dat[,4:33])
		snpdat <- data.frame("counts"=dat,"Generation"=rep(c(0,6,10),10))
		snpdat$value <- sapply(snpdat$counts, af_sync)
		anc <- filter(snpdat, Generation==0)$value
		gen6 <- filter(snpdat, Generation==6)$value - anc
		gen10 <- filter(snpdat, Generation==10)$value - anc
	  minafc <- afc
		##gen6
		if (gen6[1]>minafc){
			Rep1_6 <- c(Rep1_6,i)
			}
		if (gen6[2]>minafc){
			Rep2_6 <- c(Rep2_6,i)
			}
		if (gen6[3]>minafc){
			Rep3_6 <- c(Rep3_6,i)
			}
		if (gen6[4]>minafc){
			Rep4_6 <- c(Rep4_6,i)
			}
		if (gen6[5]>minafc){
			Rep5_6 <- c(Rep5_6,i)
			}
		if (gen6[6]>minafc){
			Rep6_6 <- c(Rep6_6,i)
			}
		if (gen6[7]>minafc){
			Rep7_6 <- c(Rep7_6,i)
			}
		if (gen6[8]>minafc){
			Rep8_6 <- c(Rep8_6,i)
			}
		if (gen6[9]>minafc){
			Rep9_6 <- c(Rep9_6,i)
			}
		if (gen6[10]>minafc){
			Rep10_6 <- c(Rep10_6,i)
			}
		##gen10
		if (gen10[1] >minafc){
			Rep1_10 <- c(Rep1_10,i)
			}
		if (gen10[2] >minafc){
			Rep2_10 <- c(Rep2_10,i)
			}
		if (gen10[3] >minafc){
			Rep3_10 <- c(Rep3_10,i)
			}
		if (gen10[4] >minafc){
			Rep4_10 <- c(Rep4_10,i)
			}
		if (gen10[5] >minafc){
			Rep5_10 <- c(Rep5_10,i)
			}
		if (gen10[6] >minafc){
			Rep6_10 <- c(Rep6_10,i)
			}
		if (gen10[7] >minafc){
			Rep7_10 <- c(Rep7_10,i)
			}
		if (gen10[8] >minafc){
			Rep8_10 <- c(Rep8_10,i)
			}
		if (gen10[9] >minafc){
			Rep9_10 <- c(Rep9_10,i)
			}
		if (gen10[10] >minafc){
			Rep10_10 <- c(Rep10_10,i)
			}
		}

	gen6_pops <- list(Rep1_6,Rep2_6,Rep3_6,Rep4_6,Rep5_6,Rep6_6,Rep7_6,Rep8_6,Rep9_6,Rep10_6)
	gen10_pops <- list(Rep1_10,Rep2_10,Rep3_10,Rep4_10,Rep5_10,Rep6_10,Rep7_10,Rep8_10,Rep9_10,Rep10_10)

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

	return(list(gen6_jac,gen10_jac))
	}

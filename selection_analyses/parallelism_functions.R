# functions to calculate Jaccard statistics and the replicate frequency spectrum (sensu Barghi et al. 2019)
# for haplotype block alleles from Stern, Anderson, Diaz, and Lee 2022
# Mostly for organizational purposes. Not great for general use

library(lme4)
library(dplyr)

jaccard_hap_blocks <- function(dat,info){
	# data = hap_blocks.rawAFC.RDS generated with prep_lmm.rawAFC.R
	# info = hab_blocks.res.RDS in data/ directory
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


		for (i in 1:length(unique(dat$loc))){
			locnum <- unique(dat$loc)[i]
			dt <- filter(dat,loc==locnum)
			#dt <- dat[which(dat$loc==snpnum)]
			gen6 <- filter(dt, Treat=="Treatment" & Generation==6)$value
			gen10 <- filter(dt, Treat=="Treatment" & Generation==10)$value
			cutoff <- info[i,]$cutoff

			##gen6
			if (gen6[1]>cutoff){
				BSE3_6 <- c(BSE3_6,dt$loc[1])
				}
			if (gen6[2]>cutoff){
				BSE4_6 <- c(BSE4_6,dt$loc[1])
				}
			if (gen6[3]>cutoff){
				BSE5_6 <- c(BSE5_6,dt$loc[1])
				}
			if (gen6[4]>cutoff){
				BSE6_6 <- c(BSE6_6,dt$loc[1])
				}
			if (gen6[5]>cutoff){
				BSE7_6 <- c(BSE7_6,dt$loc[1])
				}
			if (gen6[6]>cutoff){
				BSE8_6 <- c(BSE8_6,dt$loc[1])
				}
			if (gen6[7]>cutoff){
				BSE9_6 <- c(BSE9_6,dt$loc[1])
				}
			if (gen6[8]>cutoff){
				BSE10_6 <- c(BSE10_6,dt$loc[1])
				}
			if (gen6[9]>cutoff){
				BSE11_6 <- c(BSE11_6,dt$loc[1])
				}
			if (gen6[10]>cutoff){
				BSE12_6 <- c(BSE12_6,dt$loc[1])
				}
			##gen10
			if (gen10[1]>cutoff){
				BSE3_10 <- c(BSE3_10,dt$loc[1])
				}
			if (gen10[2]>cutoff){
				BSE4_10 <- c(BSE4_10,dt$loc[1])
				}
			if (gen10[3]>cutoff){
				BSE5_10 <- c(BSE5_10,dt$loc[1])
				}
			if (gen10[4]>cutoff){
				BSE6_10 <- c(BSE6_10,dt$loc[1])
				}
			if (gen10[5]>cutoff){
				BSE8_10 <- c(BSE8_10,dt$loc[1])
				}
			if (gen10[6]>cutoff){
				BSE9_10 <- c(BSE9_10,dt$loc[1])
				}
			if (gen10[7]>cutoff){
				BSE11_10 <- c(BSE11_10,dt$loc[1])
				}
			if (gen10[8]>cutoff){
				BSE12_10 <- c(BSE12_10,dt$loc[1])
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
						gen10_jac[pop1,pop2] <- jac
						}
					}
				}

		return(list(gen6_jac[lower.tri(gen6_jac, diag = FALSE)],gen10_jac[lower.tri(gen10_jac, diag = FALSE)]))
	}

calc_num_sel_loci_hb <- function(dat,info){
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


		for (i in 1:length(unique(dat$loc))){
			locnum <- unique(dat$loc)[i]
			dt <- filter(dat,loc==locnum)
			gen6 <- filter(dt, Treat=="Treatment" & Generation==6)$value
			gen10 <- filter(dt, Treat=="Treatment" & Generation==10)$value
			cutoff <- info[i,]$cutoff

			##gen6
			if (gen6[1]>cutoff){
				BSE3_6 <- c(BSE3_6,dt$loc[1])
				}
			if (gen6[2]>cutoff){
				BSE4_6 <- c(BSE4_6,dt$loc[1])
				}
			if (gen6[3]>cutoff){
				BSE5_6 <- c(BSE5_6,dt$loc[1])
				}
			if (gen6[4]>cutoff){
				BSE6_6 <- c(BSE6_6,dt$loc[1])
				}
			if (gen6[5]>cutoff){
				BSE7_6 <- c(BSE7_6,dt$loc[1])
				}
			if (gen6[6]>cutoff){
				BSE8_6 <- c(BSE8_6,dt$loc[1])
				}
			if (gen6[7]>cutoff){
				BSE9_6 <- c(BSE9_6,dt$loc[1])
				}
			if (gen6[8]>cutoff){
				BSE10_6 <- c(BSE10_6,dt$loc[1])
				}
			if (gen6[9]>cutoff){
				BSE11_6 <- c(BSE11_6,dt$loc[1])
				}
			if (gen6[10]>cutoff){
				BSE12_6 <- c(BSE12_6,dt$loc[1])
				}
			##gen10
			if (gen10[1]>cutoff){
				BSE3_10 <- c(BSE3_10,dt$loc[1])
				}
			if (gen10[2]>cutoff){
				BSE4_10 <- c(BSE4_10,dt$loc[1])
				}
			if (gen10[3]>cutoff){
				BSE5_10 <- c(BSE5_10,dt$loc[1])
				}
			if (gen10[4]>cutoff){
				BSE6_10 <- c(BSE6_10,dt$loc[1])
				}
			if (gen10[5]>cutoff){
				BSE8_10 <- c(BSE8_10,dt$loc[1])
				}
			if (gen10[6]>cutoff){
				BSE9_10 <- c(BSE9_10,dt$loc[1])
				}
			if (gen10[7]>cutoff){
				BSE11_10 <- c(BSE11_10,dt$loc[1])
				}
			if (gen10[8]>cutoff){
				BSE12_10 <- c(BSE12_10,dt$loc[1])
				}
			}

		gen6_pops <- list(BSE3_6,BSE4_6,BSE5_6,BSE6_6,BSE7_6,BSE8_6,BSE9_6,BSE10_6,BSE11_6,BSE12_6)
		gen10_pops <- list(BSE3_10,BSE4_10,BSE5_10,BSE6_10,BSE8_10,BSE9_10,BSE11_10,BSE12_10)

		gen6_nloci <- lapply(gen6_pops, length)
		gen10_nloci <- lapply(gen10_pops, length)
		return(list(gen6_nloci,gen10_nloci))
	}


rfs_hap_blocks <- function(dat, info){
	prop_pop6 <- c()
	prop_pop10 <- c()

	for (i in 1:length(unique(dat$loc))){
		locnum <- sort(unique(dat$loc))[i]
		dt <- filter(dat,loc==locnum)
		gen6 <- filter(dt, Treat=="Treatment" & Generation==6)$value
		gen10 <- filter(dt, Treat=="Treatment" & Generation==10)$value
		cutoff <- info[i,]$cutoff

		##gen6_sel
		gen6_sel <- gen6[gen6 >=cutoff]
		prop_pop6 <- c(prop_pop6,length(gen6_sel)/10)

		##gen10
		gen10_sel <- gen10[gen10 >=cutoff]
		prop_pop10 <- c(prop_pop10,length(gen10_sel)/8)
	}
	return(list(prop_pop6,prop_pop10))
}


af_sync <- function(counts){
	counts <- strsplit(counts,":")
	a <- as.numeric(counts[[1]][1])
	t <- as.numeric(counts[[1]][2])
	af <- a / (a+t)
	return(af)
}

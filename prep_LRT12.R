require(data.table)

# ----- Input parameters
rep <- 1:8  #replicates ids
nb_repl <- length(rep)  #number replicates
g <- c(0,4) #generations
nb_gen <- length(g)  #number generations
gen <- rep(g, nb_repl)
repl <- rep(rep, each = nb_gen)
idx_F0 <- seq(1,16,2)
idx_F4 <- seq(2,16,2)

freqs <- read.table("genobaypass_freq_noctrl",h=F)
covs <- read.table("genobaypass_cov_noctrl",h=F)
snpdet <- read.table("snpdet",h=F)

nb <- nrow(freqs)
data <- data.frame(chr = snpdet[,1], code = paste(snpdet[,2], snpdet[,3], snpdet[,4], sep = "_"), 
                     F0.BSE3 = rep("F0.BSE3", nb), cov.F0.BSE3 = covs[,idx_F0[1]], freq.F0.BSE3 = freqs[,idx_F0[1]],
                     F0.BSE4 = rep("F0.BSE4", nb), cov.F0.BSE4 = covs[,idx_F0[2]], freq.F0.BSE4 = freqs[,idx_F0[2]],
                     F0.BSE5 = rep("F0.BSE5", nb), cov.F0.BSE5 = covs[,idx_F0[3]], freq.F0.BSE5 = freqs[,idx_F0[3]],
                     F0.BSE6 = rep("F0.BSE6", nb), cov.F0.BSE6 = covs[,idx_F0[4]], freq.F0.BSE6 = freqs[,idx_F0[4]],
                     F0.BSE8 = rep("F0.BSE8", nb), cov.F0.BSE8 = covs[,idx_F0[5]], freq.F0.BSE8 = freqs[,idx_F0[5]],
                     F0.BSE9 = rep("F0.BSE9", nb), cov.F0.BSE9 = covs[,idx_F0[6]], freq.F0.BSE9 = freqs[,idx_F0[6]],
                     F0.BSE11 = rep("F0.BSE11", nb), cov.F0.BSE11 = covs[,idx_F0[7]], freq.F0.BSE11 = freqs[,idx_F0[7]],
                     F0.BSE12 = rep("F0.BSE12", nb), cov.F0.BSE12 = covs[,idx_F0[8]], freq.F0.BSE12 = freqs[,idx_F0[8]],
                     F4.BSE3 = rep("F4.BSE3", nb), cov.F4.BSE3 = covs[,idx_F4[1]], freq.F4.BSE3 = freqs[,idx_F4[1]],
                     F4.BSE4 = rep("F4.BSE4", nb), cov.F4.BSE4 = covs[,idx_F4[2]], freq.F4.BSE4 = freqs[,idx_F4[2]],
                     F4.BSE5 = rep("F4.BSE5", nb), cov.F4.BSE5 = covs[,idx_F4[3]], freq.F4.BSE5 = freqs[,idx_F4[3]],
                     F4.BSE6 = rep("F4.BSE6", nb), cov.F4.BSE6 = covs[,idx_F4[4]], freq.F4.BSE6 = freqs[,idx_F4[4]],
                     F4.BSE8 = rep("F4.BSE8", nb), cov.F4.BSE8 = covs[,idx_F4[5]], freq.F4.BSE8 = freqs[,idx_F4[5]],
                     F4.BSE9 = rep("F4.BSE9", nb), cov.F4.BSE9 = covs[,idx_F4[6]], freq.F4.BSE9 = freqs[,idx_F4[6]],
                     F4.BSE11 = rep("F4.BSE11", nb), cov.F4.BSE11 = covs[,idx_F4[7]], freq.F4.BSE11 = freqs[,idx_F4[7]],
                     F4.BSE12 = rep("F4.BSE12", nb), cov.F4.BSE12 = covs[,idx_F4[8]], freq.F4.BSE12 = freqs[,idx_F4[8]] 
                     )
  
data <- cbind(data, data.frame(x.F0.BSE3 = 2 * asin(sqrt(data$freq.F0.BSE3)), 
                                 x.F0.BSE4 = 2 * asin(sqrt(data$freq.F0.BSE4)), x.F0.BSE5 = 2 * asin(sqrt(data$freq.F0.BSE5)), 
                                 x.F0.BSE6 = 2 * asin(sqrt(data$freq.F0.BSE6)), x.F0.BSE8 = 2 * asin(sqrt(data$freq.F0.BSE8)), 
                                 x.F0.BSE9 = 2 * asin(sqrt(data$freq.F0.BSE9)), x.F0.BSE11 = 2 * asin(sqrt(data$freq.F0.BSE11)), 
                                 x.F0.BSE12 = 2 * asin(sqrt(data$freq.F0.BSE12)), x.F4.BSE3 = 2 * asin(sqrt(data$freq.F4.BSE3)), 
                                 x.F4.BSE4 = 2 * asin(sqrt(data$freq.F4.BSE4)), x.F4.BSE5 = 2 * asin(sqrt(data$freq.F4.BSE5)), 
                                 x.F4.BSE6 = 2 * asin(sqrt(data$freq.F4.BSE6)), x.F4.BSE8 = 2 * asin(sqrt(data$freq.F4.BSE8)), 
                                 x.F4.BSE9 = 2 * asin(sqrt(data$freq.F4.BSE9)), x.F4.BSE11 = 2 * asin(sqrt(data$freq.F4.BSE11)), 
                                 x.F4.BSE12 = 2 * asin(sqrt(data$freq.F4.BSE12))))
write.table(data, "lab.T1_T2_noctrl.forLRT12.z.stats.all.txt",
              sep = "\t", col.names = F, row.names = F, quote = F)


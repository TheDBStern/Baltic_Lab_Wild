require(data.table)
library(matrixStats)

# ----- Input parameters
rep <- 1:8  #replicates ids
nb_repl <- length(rep)  #number replicates
g <- c(0,10) #generations
nb_gen <- length(g)  #number generations
gen <- rep(g, nb_repl)
repl <- rep(rep, each = nb_gen)

freqs <- read.table("../lab.genobaypass_freq",h=F)
covs <- read.table("../lab.genobaypass_cov",h=F)
snpdet <- read.table("../lab.snpdet",h=F)

nb <- nrow(freqs)
data <- data.frame(chr = snpdet[,1], code = paste(snpdet[,2], snpdet[,3], snpdet[,4], sep = "_"), 
                     F0.BSE3 = rep("F0.BSE3", nb), cov.F0.BSE3 = rowSums(covs[,1:2]), freq.F0.BSE3 = rowMeans(freqs[,1:2]),
                     F0.BSE4 = rep("F0.BSE4", nb), cov.F0.BSE4 = rowSums(covs[,1:2]), freq.F0.BSE4 = rowMeans(freqs[,1:2]),
                     F0.BSE5 = rep("F0.BSE5", nb), cov.F0.BSE5 = rowSums(covs[,1:2]), freq.F0.BSE5 = rowMeans(freqs[,1:2]),
                     F0.BSE6 = rep("F0.BSE6", nb), cov.F0.BSE6 = rowSums(covs[,1:2]), freq.F0.BSE6 = rowMeans(freqs[,1:2]),
                     F0.BSE8 = rep("F0.BSE8", nb), cov.F0.BSE8 = rowSums(covs[,1:2]), freq.F0.BSE8 = rowMeans(freqs[,1:2]),
                     F0.BSE9 = rep("F0.BSE9", nb), cov.F0.BSE9 = rowSums(covs[,1:2]), freq.F0.BSE9 = rowMeans(freqs[,1:2]),
                     F0.BSE11 = rep("F0.BSE11", nb), cov.F0.BSE11 = rowSums(covs[,1:2]), freq.F0.BSE11 = rowMeans(freqs[,1:2]),
                     F0.BSE12 = rep("F0.BSE12", nb), cov.F0.BSE12 = rowSums(covs[,1:2]), freq.F0.BSE12 = rowMeans(freqs[,1:2]),
                     F10.BSE3 = rep("F10.BSE3", nb), cov.F10.BSE3 = covs[,12], freq.F10.BSE3 = freqs[,12],
                     F10.BSE4 = rep("F10.BSE4", nb), cov.F10.BSE4 = covs[,14], freq.F10.BSE4 = freqs[,14],
                     F10.BSE5 = rep("F10.BSE5", nb), cov.F10.BSE5 = covs[,16], freq.F10.BSE5 = freqs[,16],
                     F10.BSE6 = rep("F10.BSE6", nb), cov.F10.BSE6 = covs[,18], freq.F10.BSE6 = freqs[,18],
                     F10.BSE8 = rep("F10.BSE8", nb), cov.F10.BSE8 = covs[,21], freq.F10.BSE8 = freqs[,21],
                     F10.BSE9 = rep("F10.BSE9", nb), cov.F10.BSE9 = covs[,23], freq.F10.BSE9 = freqs[,23],
                     F10.BSE11 = rep("F10.BSE11", nb), cov.F10.BSE11 = covs[,26], freq.F10.BSE11 = freqs[,26],
                     F10.BSE12 = rep("F10.BSE12", nb), cov.F10.BSE12 = covs[,28], freq.F10.BSE12 = freqs[,28] 
                     )
  
data <- cbind(data, data.frame(x.F0.BSE3 = 2 * asin(sqrt(data$freq.F0.BSE3)), 
                                 x.F0.BSE4 = 2 * asin(sqrt(data$freq.F0.BSE4)), x.F0.BSE5 = 2 * asin(sqrt(data$freq.F0.BSE5)), 
                                 x.F0.BSE6 = 2 * asin(sqrt(data$freq.F0.BSE6)), x.F0.BSE8 = 2 * asin(sqrt(data$freq.F0.BSE8)), 
                                 x.F0.BSE9 = 2 * asin(sqrt(data$freq.F0.BSE9)), x.F0.BSE11 = 2 * asin(sqrt(data$freq.F0.BSE11)), 
                                 x.F0.BSE12 = 2 * asin(sqrt(data$freq.F0.BSE12)), x.F10.BSE3 = 2 * asin(sqrt(data$freq.F10.BSE3)), 
                                 x.F10.BSE4 = 2 * asin(sqrt(data$freq.F10.BSE4)), x.F10.BSE5 = 2 * asin(sqrt(data$freq.F10.BSE5)), 
                                 x.F10.BSE6 = 2 * asin(sqrt(data$freq.F10.BSE6)), x.F10.BSE8 = 2 * asin(sqrt(data$freq.F10.BSE8)), 
                                 x.F10.BSE9 = 2 * asin(sqrt(data$freq.F10.BSE9)), x.F10.BSE11 = 2 * asin(sqrt(data$freq.F10.BSE11)), 
                                 x.F10.BSE12 = 2 * asin(sqrt(data$freq.F10.BSE12))))
write.table(data, "lab.T0_T2.sel.z.stats.all.txt",
              sep = "\t", col.names = F, row.names = F, quote = F)


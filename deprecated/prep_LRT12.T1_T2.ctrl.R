require(data.table)
library(matrixStats)

# ----- Input parameters
rep <- 1:2  #replicates ids
nb_repl <- length(rep)  #number replicates
g <- c(4,10) #generations
nb_gen <- length(g)  #number generations
gen <- rep(g, nb_repl)
repl <- rep(rep, each = nb_gen)

freqs <- read.table("../lab.genobaypass_freq",h=F)
covs <- read.table("../lab.genobaypass_cov",h=F)
snpdet <- read.table("../lab.snpdet",h=F)

nb <- nrow(freqs)
data <- data.frame(chr = snpdet[,1], code = paste(snpdet[,2], snpdet[,3], snpdet[,4], sep = "_"), 
                     F6.BSE1 = rep("F6.BSE1", nb), cov.F6.BSE1 = covs[,5], freq.F6.BSE1 = freqs[,5],
                     F6.BSE2 = rep("F6.BSE2", nb), cov.F6.BSE2 = covs[,8], freq.F6.BSE2 = freqs[,8],
                     F10.BSE1 = rep("F10.BSE1", nb), cov.F10.BSE1 = rowSums(covs[,6:7]), freq.F10.BSE1 = rowMeans(freqs[,6:7]),
                     F10.BSE2 = rep("F10.BSE2", nb), cov.F10.BSE2 = rowSums(covs[,9:10]), freq.F10.BSE2 = rowMeans(freqs[,9:10])
                     )
  
data <- cbind(data, data.frame(x.F6.BSE1 = 2 * asin(sqrt(data$freq.F6.BSE1)), 
                                 x.F6.BSE2 = 2 * asin(sqrt(data$freq.F6.BSE2)),
                                 x.F10.BSE1 = 2 * asin(sqrt(data$freq.F10.BSE1)), 
                                 x.F10.BSE2 = 2 * asin(sqrt(data$freq.F10.BSE2))))
write.table(data, "lab.T1_T2.ctrl.z.stats.all.txt",
              sep = "\t", col.names = F, row.names = F, quote = F)


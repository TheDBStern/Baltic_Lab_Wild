library(dplyr)

dat <- read.table(gzfile("../wild.noPBE.genobaypass_treemix.gz"),h=T)

for (i in 1:100){
	dat_s <- sample_n(dat,751)
	write.table(dat_s, file=paste("tmp",i,sep=''),quote=F,sep=' ',row.names=F,)
	system(paste("gzip tmp",i,sep=''))
	
	system(paste('threepop -k 1 -i tmp',i,'.gz > threepop_rep',i,sep=''))
	system(paste('fourpop -k 1 -i tmp',i,'.gz > fourpop_rep',i,sep=''))
}


##### summarize results

###F3
dat <- read.table('threepop_results.clean.txt',h=F)
colnames(dat) <- c("Pops","F3","StdErr","Zscore")

#full summary
F3_sum <- dat %>%
			group_by(Pops) %>%
			summarize(avg = mean(F3)) %>%
			arrange(-avg)
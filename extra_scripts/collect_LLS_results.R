library(matrixStats)

snpdet <- read.table('../lab.snpdet',h=F)

BSE1 <- read.table('BSE1.PplanII.txt',h=T)
BSE2 <- read.table('BSE2.PplanII.txt',h=T)
BSE3 <- read.table('BSE3.PplanII.txt',h=T)
BSE4 <- read.table('BSE4.PplanII.txt',h=T)
BSE5 <- read.table('BSE5.PplanII.txt',h=T)
BSE6 <- read.table('BSE6.PplanII.txt',h=T)
BSE8 <- read.table('BSE8.PplanII.txt',h=T)
BSE9 <- read.table('BSE9.PplanII.txt',h=T)
BSE11 <- read.table('BSE11.PplanII.txt',h=T)
BSE12 <- read.table('BSE12.PplanII.txt',h=T)

sel <- as.matrix(data.frame('BSE3'=BSE3$s,'BSE4'=BSE4$s,'BSE5'=BSE5$s,'BSE6'=BSE6$s,'BSE8'=BSE8$s,'BSE9'=BSE9$s,'BSE11'=BSE11$s,'BSE12'=BSE12$s))
sel <- cbind(sel,"Mean"=rowMeans(sel,na.rm=T))
sel <- cbind(sel,"Min"=rowMins(sel[,1:8],na.rm=T))
sel <- cbind(sel,"Max" =rowMaxs(sel[,1:8],na.rm=T))
sel <- cbind(sel,"Median" =rowMedians(sel[,1:8],na.rm=T))

sel <- data.frame('Pseudo_Transcript'=snpdet[,1],'Pseudo_Position'=snpdet[,2],sel$Mean)

sig_pos <- sel[which(sel$Mean > 0 & sel$Min > 0 & sel$Max > 0),]
sig_neg <- sel[which(sel$Mean < 0 & sel$Min < 0 & sel$Max < 0),]

sel_sig <- rbind(sig_pos,sig_neg)




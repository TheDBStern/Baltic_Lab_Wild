dat <- read.table('lab.genobaypass_freqs',h=F)
snpdet <- read.table('lab.snpdet',h=F)

dat$Transcript <- snpdet[,1]
dat$Position <- snpdet[,2]

#filter any row with a population with 0 AF
row_sub = apply(dat[,1:28], 1, function(row) all(row !=0 ))
filtered <- dat[row_sub,]

#filter any row with a population with 1 AF
row_sub = apply(filtered[,1:28], 1, function(row) all(row !=1 ))
filtered <- filtered[row_sub,]

write.table(filtered, 'lab.genobaypass_freqs_filtered',quote=F, row.names=F, sep=' ')
write.table(filtered[,29:30], 'lab.snpdet_filtered',quote=F, row.names=F, sep=' ')

## filter the genobaypass file

geno <- read.table('lab.genobaypass',h=F)

geno$Transcript <- snpdet[,1]
geno$Position <- snpdet[,2]

geno_filt <- merge(geno,filtered, by=c("Transcript","Position"))
write.table(geno_filt[,2:58], 'lab.genobaypass_filtered',quote=F, row.names=F, sep=' ')


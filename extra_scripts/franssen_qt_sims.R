library(dplyr)

sig <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.new.cmh05.lrt05.RDS')
#sig <- sample_n(sig_all,25)

starting_afs <- paste(round(sig$T0_AF_rising,3),collapse=',')
#starting_afs <- paste(rep(0.05,10),collapse=',')
#sum of effects = 137.186
#contrs <- round(abs(sig$estS_lm),3)
#test_contrs <- rep(4,50)
#contrs_w_test <- paste(c(test_contrs,round(abs(sig$estS_lm[51:1156]),3)),collapse=',')
contrs <- paste(round(abs(sig$selCoef),3),collapse=',')
#contrs <- paste(rep(1,10),collapse=',')
sel <- '0.5:6.5:0.9:0.2'
pre <- 'opt0.9.sd0.3.ne2000.1156loci'
gens <- '6,10'

cmd <- paste('python ~/Desktop/Genomics_Programs/Franssen_qt_sim/python_qt_simulator/frequencyAt-pheno-quantitative.py --Ne 2000 --selection ',
			sel,' --variant-contributions ',contrs,' -p ',starting_afs,' --snapshots ',gens,' --repeat-simulations 10 --outfile-freq ',
			pre,'.freqs.out --outfile-pheno ',pre,'.pheno.out',sep='')

system(cmd)

## testing number of loci
for (loc in 1156){
		for (rep in 1:10){
			#sig <- sample_n(sig_all,loc)
			starting_afs <- paste(round(sig$T0_AF_rising,3),collapse=',')
			contrs <- paste(round(abs(sig$selCoef),3),collapse=',')
			sel <- '0.5:4.5:0.9:0.3'
			pre <- paste(c('opt0.9.sd0.3.ne2000.loc',loc,'.rep',rep),collapse='')
			print(pre)
			gens <- '6,10'
			cmd <- paste('python ~/Desktop/Genomics_Programs/Franssen_qt_sim/python_qt_simulator/frequencyAt-pheno-quantitative.py --Ne 2000 --selection ',
				sel,' --variant-contributions ',contrs,' -p ',starting_afs,' --snapshots ',gens,' --repeat-simulations 10 --outfile-freq ',
				pre,'.freqs.out --outfile-pheno ',pre,'.pheno.out',sep='')
			system(cmd)
			}
		}

#barghi 0.5:4.5:0.6:0.3
#franssen 0.5:1.5:0.4:0.2
#test1 0.5:4.5:0.95:0.1
#test2

## calculate jaccard
library(dplyr)
library(stringr)
simdat <- read.table('opt0.9.sd0.3.ne2000.loc40.rep4.freqs.out',h=F)
colnames(simdat) <- c("snp_id","generation","frequency","replicate")

res <- calc_jaccard(simdat)

## loop through all files, calculate jaccard, put into table
plotdat <- data.frame("Loci"=c(),"Generation"=c(),"Jaccard"=c())


files <- list.files(path = ".", pattern = 'freqs')
for (i in 1:length(files)){
	print(i)
	file <- files[i]
	loci <- as.numeric(str_remove(strsplit(file,'\\.')[[1]][6],"loc"))
	simdat <- read.table(file,h=F)
	colnames(simdat) <- c("snp_id","generation","frequency","replicate")
	res <- calc_jaccard(simdat)
	df <- data.frame("Loci"=rep(loci,180),"Generation"=c(rep("Asix",90),rep("Ten",90)),"Jaccard"=c(res[[1]],res[[2]]))
	plotdat <- rbind(plotdat,df)
	}

p <- ggplot(plotdat, aes(x = Loci, y = Jaccard, color=Generation,fill=Generation)) 
p +   stat_summary(fun=mean, geom="point", position=position_dodge(width=1)) +
	stat_smooth(method="loess",se = TRUE) +
	scale_color_manual(values=c("#707070","#c4c4c4")) + 
	scale_fill_manual(values=c("#707070","#c4c4c4")) +
	theme_classic()



calc_jaccard <- function(simdat){
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
  	minafc <- 0.01
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
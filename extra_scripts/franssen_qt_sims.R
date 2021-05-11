library(dplyr)

sig_all <- readRDS('lab.sig.RDS')
all <- readRDS('lab.all.RDS')


## testing number of loci
for (loc in seq(10,200,10)){
		for (rep in 1:100){
			sig <- sample_n(all,loc)
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

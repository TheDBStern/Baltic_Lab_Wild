## polygenic adaptation figure code

library(ggplot2)
library(dplyr)

sal <- c(2.5,2.6,3.2,6,6,0,5.63,0,5.51,5.06,7,6.1)
bin <- c(1,1,1,0,0,1,0,1,0,0,0,0)
lats <- c(65.0167,65.3836,63.5169,59.8822,59.8822,52.699,54.3333,52.574,57.5844,57.4342,51.302,58.8553)
longs <- c(22.6833,23.4706,21.5217,23.2539,23.2539,5.29,10.15,5.033,22.9439,23.7097,4.286,17.6283)
pop_names <- c("BB1E","BB2E","GBE","HF1E","HF2E","IJE","KIE","MME","RG1E","RG2E","SCE","STE")
pop_names <- c("BB1E","BB2E","GBE","HF1E","HF2E","KIE","RG1E","RG2E","STE")

region <- c(rep("Baltic",5),"North_Sea","Baltic","North_Sea","Baltic","Baltic","North_Sea","Baltic")
gendist <- c(0.617,0.613,0.608,0.612,0.615,0.187,1,-0.134,0.578,0.608,-0.161,0.607)

load('all_sigsnps_baltic/Output/genetic.values.Robj') #gvs
load('all_sigsnps_baltic/Output/asymptotic.pVals.Robj') #asymp.p.vals
load('all_sigsnps_baltic/Output/nullStats.Robj') #null.stats
load('all_sigsnps_baltic/Output/theStats.Robj') #the.stats
load('all_sigsnps_baltic/Output/pVals.Robj') #p.vals

gv_dat <- data.frame(Gen.val =gvs,sal,bin,lats,longs,pop_names,region,gendist)
gv_dat <- gv_dat[order(gv_dat$region),]

#plot genetic values by env. variable
p <- ggplot(gv_dat,aes(x=sal,y=Gen.val, fill=region,color=region))
p + geom_point(size=8) + 
	scale_color_manual(values=c("Baltic"="#ffc469","North_Sea"="#7501BC")) +
	theme_bw()

qx_null <- data.frame("qx_null"=null.stats$Qx)

qx <- data.frame("Qx"=the.stats$Qx)

p <- ggplot(data=qx_null,aes(x = qx_null))
p + geom_histogram(bins=50)+
	geom_vline(data=qx, aes(xintercept=Qx,  colour="red"),
               linetype="dashed", size=1) +
	theme_classic()

bardat <- data.frame("component"=c("Qx","LD","Fst"),"value"=c(the.stats$Qx,the.stats$LD.component,the.stats$Fst.comp))

ggplot(data=bardat, aes(fill=component, y=value, x=component)) + 
    geom_bar(position="dodge", stat="identity") + theme_bw()
    

## baltic only

library(ggplot2)
library(dplyr)

sal <- c(2.5,2.6,3.2,6,6,5.63,5.51,5.06,6.1)
bin <- c(1,1,1,0,0,0,0,0,0)
lats <- c(65.0167,65.3836,63.5169,59.8822,59.8822,54.3333,57.5844,57.4342,58.8553)
longs <- c(22.6833,23.4706,21.5217,23.2539,23.2539,10.15,22.9439,23.7097,17.6283)
pop_names <- c("BB1E","BB2E","GBE","HF1E","HF2E","KIE","RG1E","RG2E","STE")
gendist <- c(0.617,0.613,0.608,0.612,0.615,1,0.578,0.608,0.607)

load('all_sigsnps_baltic_only/Output/genetic.values.Robj') #gvs
load('all_sigsnps_baltic_only/Output/asymptotic.pVals.Robj') #asymp.p.vals
load('all_sigsnps_baltic_only/Output/nullStats.Robj') #null.stats
load('all_sigsnps_baltic_only/Output/theStats.Robj') #the.stats
load('all_sigsnps_baltic_only/Output/pVals.Robj') #p.vals

gv_dat <- data.frame(Gen.val =gvs,sal,bin,lats,longs,pop_names,gendist)

#plot genetic values by env. variable
p <- ggplot(gv_dat,aes(x=sal,y=Gen.val))
p + geom_point(size=8) + 
	theme_bw()

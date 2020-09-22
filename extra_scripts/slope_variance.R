library(data.table)
library(dplyr)

dat <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/full_data.rawFreqs.RDS')

#bse0_tmp <- filter(dat,Beaker=="BSE-0")
#bse0_tmp$Treat <- rep("Treatment",nrow(bse0_tmp))
#bse0b_tmp <- filter(dat,Beaker=="BSE-0B")
#bse0b_tmp$Treat <- rep("Treatment",nrow(bse0b_tmp))
#dat <- rbind(dat,bse0_tmp)
#dat <- rbind(dat,bse0b_tmp)
dat[which(dat$value==0),3] <- 0.01
dat[which(dat$value==1),3] <- 0.99

dat$AF_logit <- log(dat$value/(1-dat$value))
dat$asin_sqrt <- 2*asin(sqrt(dat$value))

sig <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/lab.res.new.cmh05.lrt05.RDS')
colnames(sig)[1:2] <- c('Transcript','Position')
sig_dat <- merge(dat,sig,by=c('Transcript','Position'))

## 
# 1
# estimate slope for each line, and get variance in slope across lines
##

lm_test <-function(i){
  print(i)
  SNPdat= dat[which(dat$SNP == i),]
  SNPdat <- filter(SNPdat,Treat=="Treatment" | variable=="Ancestor")
  anc <- data.frame(SNP=i,value=mean(filter(SNPdat,variable=="Ancestor")$value), Generation=0,Coverage=sum(filter(SNPdat,variable=="Ancestor")$Coverage))
  anc$AF_logit <- log(anc$value/(1-anc$value))
  anc$asin_sqrt <- 2*asin(sqrt(anc$value))
  beakers <- c("BSE-3","BSE-4","BSE-5","BSE-6","BSE-8","BSE-9","BSE-11","BSE-12")
  slope <- c()
  for (x in 1:length(beakers)){
  	beak <- beakers[[x]]
  	beakdat <- filter(SNPdat,Beaker==beak)
  	snpdat <- rbind(anc,beakdat[c(1,3,4,6,10,11)])
  	neff <- (100*snpdat$Coverage-1)/ (100+snpdat$Coverage)
    mod <- lm(AF_logit~Generation, weights=neff, data=snpdat)
    sum <- summary(mod)
    slope <- c(slope,coef(sum)[2,1])
  }
  return(var(slope))
}

lm_res <- mclapply(1:length(sig_snp_nums),lm_test,mc.cores=4)
lm_res <- lapply(1:length(unique(dat$SNP)),lm_test)

lm_res <- unlist(lm_res)


##
# 2
# estimate snp frequency variance and covariance for each line
##

cov_test <-function(i){
  print(i)
  SNPdat= dat[which(dat$SNP == i),]
  SNPdat <- filter(SNPdat,Treat=="Treatment" | variable=="Ancestor")
  beakers <- c("BSE-3","BSE-4","BSE-5","BSE-6","BSE-8","BSE-9","BSE-11","BSE-12")
  freqdat <- data.frame()
  for (x in 1:length(beakers)){
  	beak <- beakers[[x]]
  	beakdat <- filter(SNPdat,Beaker==beak)
  	df <- data.frame(F6=filter(beakdat,Generation==6)$asin_sqrt, F10=filter(beakdat,Generation==10)$asin_sqrt)
  	freqdat <- rbind(freqdat,df)
  }  
  t.freqdat <- t(freqdat)
  covres <- cov(t.freqdat)
  return(mean(covres[lower.tri(covres)]))
}

cov_res <- mclapply(1:length(unique(dat$SNP)),cov_test,mc.cores=3)
cov_res <- unlist(cov_res)

##
# 3
# is snp freq variance among treatment lines correlated with selection coefficient?
##
sig_dat_gen10 <- filter(sig_dat,Treat=="Treatment" & Generation=="10")

var_res <- sig_dat_gen10 %>%
		group_by(Transcript,Position) %>%
		summarize(var=sqrt(var(AF_logit)))

####
plotdat <- data.frame("cov"=c(sig$mean_cov,res_all$mean_cov),"Group"=c(rep("Selected",nrow(sig)),rep("Non-Selected",nrow(res_all))))

p <- ggplot(plotdat, aes(x = Group, y=cov, color=Group,fill=Group))
p + geom_boxplot()

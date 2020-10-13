## power analysis for new selection experiment
#d: Effect size (Cohen's d) - difference between the means
          # divided by the pooled standard deviation
#sig.level: Significance level (Type I error probability)
#power: Power of test (1 minus Type II error probability)

library(dplyr)
library(pwr)
library(ssize)
#load baltic sea data
dat <- readRDS('/Users/dbstern/Desktop/Baltic_sea_project/pseudoref/varscan/full_data.rawFreqs.RDS')
dat$angular <- 1*asin(sqrt(dat$value))
anc <- filter(dat,Generation==0)
treat_gen10 <- filter(dat,Generation==10 & Treat=="Treatment")

gen10 <- filter(dat,Generation==10)
start_end <- rbind(anc, treat_gen10)

# calculate standard deviation in SNP frequency at gen10 for raw freqs
sd_res <- start_end %>%
      group_by(SNP) %>%
      summarize(sd_raw= sd(value), sd_ang = sd(angular))


pwrt<-pwr.t.test(d=0.01/mean(sd_res$sd_raw),n=c(5,6,7,8,9,10),sig.level=.05,type="two.sample",alternative="two.sided")

plot(pwrt$n,pwrt$power,type="b",xlab="sample size",ylab="power")

#
pwrt.fdr <- power.t.test.FDR(sd=mean(sd_res$sd_raw), n=5, delta=NULL,
                 FDR.level=0.05,
                 pi0=0.99,
                 power=0.8,
                 type="two.sample",
                 alternative="two.sided" )

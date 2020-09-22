library(data.table)
library(dplyr)

# simulated data had 10 selected lines, 3 timepoints each
# format start, gen6, gen10 for 10 selected lines


# read in data
genoF <- as.data.table(read.table("emp_vals.1156sel.selStart.rep10.freq", h=T, stringsAsFactors=F))
genoF$SNP<- 1:nrow(genoF)
genoF_cov <- as.data.table(read.table("emp_vals.1156sel.selStart.rep10.cov", h=T, stringsAsFactors=F))
#snpdet <- read.table("lab.snpdet",h=F)[,1:2]

## reformat
setDT(genoF)
MgenoF<-melt(genoF, id.vars ="SNP")
MgenoF$variable<-as.character(MgenoF$variable)
MgenoF$Generation<-rep(NA,nrow(MgenoF))
MgenoF$Treat<-rep(NA,nrow(MgenoF))
MgenoF<- MgenoF %>%
  mutate(variable = replace(variable, variable=="V1", "Treatment0")) %>%
  mutate(variable = replace(variable, variable=="V2", "Treatment6")) %>%
  mutate(variable = replace(variable, variable=="V3", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V4", "Treatment0")) %>% 
  mutate(variable = replace(variable, variable=="V5", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V6", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V7", "Treatment0")) %>% 
  mutate(variable = replace(variable, variable=="V8", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V9", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V10", "Treatment0")) %>% 
  mutate(variable = replace(variable, variable=="V11", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V12", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V13", "Treatment0")) %>% 
  mutate(variable = replace(variable, variable=="V14", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V15", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V16", "Treatment0")) %>% 
  mutate(variable = replace(variable, variable=="V17", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V18", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V19", "Treatment0")) %>% 
  mutate(variable = replace(variable, variable=="V20", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V21", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V22", "Treatment0")) %>% 
  mutate(variable = replace(variable, variable=="V23", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V24", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V25", "Treatment0")) %>% 
  mutate(variable = replace(variable, variable=="V26", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V27", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V28", "Treatment0")) %>% 
  mutate(variable = replace(variable, variable=="V29", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V30", "Treatment10")) %>%  
  mutate(Treat = replace(Treat, variable=="Treatment0", "Treatment"),  Generation  = replace(Generation, variable=="Treatment0", "0"))  %>% 
  mutate(Treat= replace(Treat, variable=="Treatment6", "Treatment"),  Generation  = replace(Generation, variable=="Treatment6", "6"))  %>% 
  mutate(Treat = replace(Treat, variable=="Treatment10", "Treatment"),  Generation  = replace(Generation, variable=="Treatment10", "10")) 

## Adding coverage of each
cov<-genoF_cov
setDT(cov)
cov<-melt(cov)
MgenoF$Coverage <- cov$value

## Adding which beaker
beaker<-c(rep("T1", 3*nrow(genoF)),rep("T2", 3*nrow(genoF)),rep("T3", 3*nrow(genoF)),rep("T4", 3*nrow(genoF)),rep("T5", 3*nrow(genoF)),rep("T6", 3*nrow(genoF)),rep("T7", 3*nrow(genoF)),rep("T8", 3*nrow(genoF)),rep("T9", 3*nrow(genoF)),rep("T10", 3*nrow(genoF)))
MgenoF$Beaker<-beaker


## adjust data types and save
MgenoF$Generation <- as.numeric(MgenoF$Generation)
MgenoF$Coverage <- as.numeric(MgenoF$Coverage)
MgenoF$Beaker <- as.factor(MgenoF$Beaker)
MgenoF$Treat <- as.factor(MgenoF$Treat)
MgenoF <- as.data.table(MgenoF)

saveRDS(MgenoF,'full_data_lmm.sim.1156sel.raw.rep10.RDS')
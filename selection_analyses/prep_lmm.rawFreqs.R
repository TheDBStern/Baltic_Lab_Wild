library(data.table)
library(dplyr)

# read in data
genoF <- as.data.table(read.table("data/hap_blocks.freq", h=F, stringsAsFactors=F))
genoF$loc<- 1:nrow(genoF)
genoF_cov <- as.data.table(read.table("data/hap_blocks.cov", h=F, stringsAsFactors=F))

## reformat
setDT(genoF)
MgenoF<-melt(genoF, id.vars ="loc")
MgenoF$variable<-as.character(MgenoF$variable)
MgenoF$Generation<-rep(NA,nrow(MgenoF))
MgenoF$Treat<-rep(NA,nrow(MgenoF))
MgenoF<- MgenoF %>%
  mutate(variable = replace(variable, variable=="V1", "Ancestor")) %>%
  mutate(variable = replace(variable, variable=="V2", "Ancestor")) %>%
  mutate(variable = replace(variable, variable=="V3", "Control20")) %>%
  mutate(variable = replace(variable, variable=="V4", "Control20")) %>%
  mutate(variable = replace(variable, variable=="V5", "Control6")) %>%
  mutate(variable = replace(variable, variable=="V6", "Control10")) %>%
  mutate(variable = replace(variable, variable=="V7", "Control10")) %>%
  mutate(variable = replace(variable, variable=="V8", "Control6")) %>%
  mutate(variable = replace(variable, variable=="V9", "Control10")) %>%
  mutate(variable = replace(variable, variable=="V10", "Control10")) %>%
  mutate(variable = replace(variable, variable=="V11", "Treatment6")) %>%
  mutate(variable = replace(variable, variable=="V12", "Treatment10")) %>%
  mutate(variable = replace(variable, variable=="V13", "Treatment6")) %>%
  mutate(variable = replace(variable, variable=="V14", "Treatment10")) %>%
  mutate(variable = replace(variable, variable=="V15", "Treatment6")) %>%
  mutate(variable = replace(variable, variable=="V16", "Treatment10")) %>%
  mutate(variable = replace(variable, variable=="V17", "Treatment6")) %>%
  mutate(variable = replace(variable, variable=="V18", "Treatment10")) %>%
  mutate(variable = replace(variable, variable=="V19", "Treatment6")) %>%
  mutate(variable = replace(variable, variable=="V20", "Treatment6")) %>%
  mutate(variable = replace(variable, variable=="V21", "Treatment10")) %>%
  mutate(variable = replace(variable, variable=="V22", "Treatment6")) %>%
  mutate(variable = replace(variable, variable=="V23", "Treatment10")) %>%
  mutate(variable = replace(variable, variable=="V24", "Treatment6")) %>%
  mutate(variable = replace(variable, variable=="V25", "Treatment6")) %>%
  mutate(variable = replace(variable, variable=="V26", "Treatment10")) %>%
  mutate(variable = replace(variable, variable=="V27", "Treatment6")) %>%
  mutate(variable = replace(variable, variable=="V28", "Treatment10")) %>%
  mutate(Treat = replace(Treat, variable=="Ancestor", "Control"),  Generation  = replace(Generation, variable=="Ancestor", "0")) %>%
  mutate(Treat = replace(Treat, variable=="Control20", "Control"),  Generation  = replace(Generation, variable=="Control20", "20")) %>%
  mutate(Treat = replace(Treat, variable=="Control6", "Control"),  Generation  = replace(Generation, variable=="Control6", "6")) %>%
  mutate(Treat = replace(Treat, variable=="Control10", "Control"),  Generation  = replace(Generation, variable=="Control10", "10"))  %>%
  mutate(Treat= replace(Treat, variable=="Treatment6", "Treatment"),  Generation  = replace(Generation, variable=="Treatment6", "6"))  %>%
  mutate(Treat = replace(Treat, variable=="Treatment10", "Treatment"),  Generation  = replace(Generation, variable=="Treatment10", "10"))

## Adding coverage of each
cov<-genoF_cov
setDT(cov)
cov<-melt(cov)
MgenoF$Coverage <- cov$value

## Adding which beaker
beaker<-c(rep("BSE-0", nrow(genoF)),rep("BSE-0B", nrow(genoF)),rep("BS3C", nrow(genoF)),rep("BS4C", nrow(genoF)),rep("BSE-1", 3*nrow(genoF)), rep("BSE-2", 3*nrow(genoF)), rep("BSE-3", 2*nrow(genoF)), rep("BSE-4", 2*nrow(genoF)),rep("BSE-5", 2*nrow(genoF)), rep("BSE-6", 2*nrow(genoF)),rep("BSE-7", nrow(genoF)),rep("BSE-8", 2*nrow(genoF)),rep("BSE-9", 2*nrow(genoF)),rep("BSE-10", nrow(genoF)),rep("BSE-11", 2*nrow(genoF)),rep("BSE-12", 2*nrow(genoF)))
MgenoF$Beaker<-beaker



## adjust data types and save
MgenoF$Generation <- as.numeric(MgenoF$Generation)
MgenoF$Coverage <- as.numeric(MgenoF$Coverage)
MgenoF$Beaker <- as.factor(MgenoF$Beaker)
MgenoF$Treat <- as.factor(MgenoF$Treat)
MgenoF <- as.data.table(MgenoF)

saveRDS(MgenoF,'hap_blocks.rawFreqs.RDS')

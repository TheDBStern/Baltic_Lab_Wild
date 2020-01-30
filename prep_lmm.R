library(data.table)
library(dplyr)

# read in data
genoF <- read.table("genobaypass_lab_freq", h=F, stringsAsFactors=F)
genoF$SNP<- 1:nrow(genoF)
genoF_cov <- read.table("genobaypass_lab_cov", h=F, stringsAsFactors=F)

## calculate divergence from ancestor using Kelly and Hughes 2019 angular transformation
data=select(genoF, V1,V2)
genoF$Ancestral= 2*asin(sqrt(rowMeans(data)))

genoF<-genoF %>% 
mutate(V3= 2*asin(sqrt(V3))-Ancestral) %>% 
mutate(V4= 2*asin(sqrt(V4))-Ancestral) %>% 
mutate(V5= 2*asin(sqrt(V5))-Ancestral)%>% 
mutate(V6= 2*asin(sqrt(V6))-Ancestral)%>% 
mutate(V7= 2*asin(sqrt(V7))-Ancestral)%>% 
mutate(V8= 2*asin(sqrt(V8))-Ancestral)%>% 
mutate(V9= 2*asin(sqrt(V9))-Ancestral)%>% 
mutate(V10= 2*asin(sqrt(V10))-Ancestral)%>% 
mutate(V11= 2*asin(sqrt(V11))-Ancestral)%>% 
mutate(V12= 2*asin(sqrt(V12))-Ancestral)%>% 
mutate(V13= 2*asin(sqrt(V13))-Ancestral)%>% 
mutate(V14= 2*asin(sqrt(V14))-Ancestral)%>% 
mutate(V15= 2*asin(sqrt(V15))-Ancestral)%>% 
mutate(V16= 2*asin(sqrt(V16))-Ancestral)%>% 
mutate(V17= 2*asin(sqrt(V17))-Ancestral)%>% 
mutate(V18= 2*asin(sqrt(V18))-Ancestral)%>% 
mutate(V19= 2*asin(sqrt(V19))-Ancestral)%>% 
mutate(V20= 2*asin(sqrt(V20))-Ancestral)%>% 
mutate(V21= 2*asin(sqrt(V21))-Ancestral)%>% 
mutate(V22= 2*asin(sqrt(V22))-Ancestral)%>% 
mutate(V23= 2*asin(sqrt(V23))-Ancestral)%>% 
mutate(V24= 2*asin(sqrt(V24))-Ancestral)%>% 
mutate(V25= 2*asin(sqrt(V25))-Ancestral)%>% 
mutate(V26= 2*asin(sqrt(V26))-Ancestral)%>% 
mutate(V27= 2*asin(sqrt(V27))-Ancestral)%>% 
mutate(V28= 2*asin(sqrt(V28))-Ancestral)

## reformat
genoF<-select(genoF,3:29)
MgenoF<-melt(genoF, id.vars ="SNP")
MgenoF$variable<-as.character(MgenoF$variable)
MgenoF$Generation<-rep(NA,nrow(MgenoF))
MgenoF$Treat<-rep(NA,nrow(MgenoF))
MgenoF<- MgenoF %>% 
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
  mutate(Treat = replace(Treat, variable=="Control20", "Control"),  Generation  = replace(Generation, variable=="Control20", "20")) %>% 
  mutate(Treat = replace(Treat, variable=="Control6", "Control"),  Generation  = replace(Generation, variable=="Control6", "6")) %>% 
  mutate(Treat = replace(Treat, variable=="Control10", "Control"),  Generation  = replace(Generation, variable=="Control10", "10"))  %>% 
  mutate(Treat= replace(Treat, variable=="Treatment6", "Treatment"),  Generation  = replace(Generation, variable=="Treatment6", "6"))  %>% 
  mutate(Treat = replace(Treat, variable=="Treatment10", "Treatment"),  Generation  = replace(Generation, variable=="Treatment10", "10")) 

## Adding coverage of each
cov<-select(genoF_cov,3:28)
cov<-melt(cov)
MgenoF$Coverage <- cov$value

## Adding which beaker
beaker<-c(rep("BS3C", nrow(genoF)),rep("BS4C", nrow(genoF)),rep("BSE-1", 3*nrow(genoF)), rep("BSE-2", 3*nrow(genoF)), rep("BSE-3", 2*nrow(genoF)), rep("BSE-4", 2*nrow(genoF)),rep("BSE-5", 2*nrow(genoF)), rep("BSE-6", 2*nrow(genoF)),rep("BSE-7", nrow(genoF)),rep("BSE-8", 2*nrow(genoF)),rep("BSE-9", 2*nrow(genoF)),rep("BSE-10", nrow(genoF)),rep("BSE-11", 2*nrow(genoF)),rep("BSE-12", 2*nrow(genoF)))
MgenoF$Beaker<-beaker

## adjust data types and save
MgenoF$Generation <- as.numeric(MgenoF$Generation)
MgenoF$Coverage <- as.numeric(MgenoF$Coverage)
MgenoF$Beaker <- as.factor(MgenoF$Beaker)
MgenoF$Treat <- as.factor(MgenoF$Treat)
MgenoF <- as.data.table(MgenoF)

saveRDS(MgenoF,'full_data_lmm.RDS')
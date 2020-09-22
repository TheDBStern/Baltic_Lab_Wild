library(data.table)
library(dplyr)

# simulated data had 100k snps with 4 control lines, 10 selected lines, 3 timepoints each

# read in data
genoF <- as.data.table(read.table("emp_vals.1ksel_1knon.empNe.selstart.freq", h=T, stringsAsFactors=F))
genoF$SNP<- 1:nrow(genoF)
genoF_cov <- as.data.table(read.table("emp_vals.1ksel_1knon.empNe.selstart.cov", h=T, stringsAsFactors=F))

## calculate divergence from ancestor using Kelly and Hughes 2019 angular transformation
anc=select(genoF, V1,V4,V7,V10,V13,V16,V19,V22,V25,V28,V31,V34,V37,V40)
genoF$Ancestral= rowMeans(anc)

genoF<-genoF %>% 
mutate(V2= V2-Ancestral) %>% 
mutate(V3= V3-Ancestral) %>% 
mutate(V5= V5-Ancestral) %>% 
mutate(V6= V6-Ancestral)%>% 
mutate(V8= V8-Ancestral)%>% 
mutate(V9= V9-Ancestral)%>% 
mutate(V11= V11-Ancestral)%>% 
mutate(V12= V12-Ancestral)%>% 
mutate(V14= V14-Ancestral)%>% 
mutate(V15= V15-Ancestral)%>% 
mutate(V17= V17-Ancestral)%>% 
mutate(V18= V18-Ancestral)%>% 
mutate(V20= V20-Ancestral)%>% 
mutate(V21= V21-Ancestral)%>% 
mutate(V23= V23-Ancestral)%>% 
mutate(V24= V24-Ancestral)%>% 
mutate(V26= V26-Ancestral)%>%
mutate(V27= V27-Ancestral)%>% 
mutate(V29= V29-Ancestral)%>% 
mutate(V30= V30-Ancestral)%>% 
mutate(V32= V32-Ancestral)%>% 
mutate(V33= V33-Ancestral)%>% 
mutate(V35= V35-Ancestral)%>% 
mutate(V36= V36-Ancestral)%>% 
mutate(V38= V38-Ancestral)%>% 
mutate(V39= V39-Ancestral)%>% 
mutate(V41= V41-Ancestral)%>% 
mutate(V42= V42-Ancestral)

## reformat
setDT(genoF)
genoF<-select(genoF,2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33,35,36,38,39,41,42,43)
MgenoF<-melt(genoF, id.vars ="SNP")
MgenoF$variable<-as.character(MgenoF$variable)
MgenoF$Generation<-rep(NA,nrow(MgenoF))
MgenoF$Treat<-rep(NA,nrow(MgenoF))
MgenoF<- MgenoF %>% 
  mutate(variable = replace(variable, variable=="V2", "Control6")) %>% 
  mutate(variable = replace(variable, variable=="V3", "Control10")) %>% 
  mutate(variable = replace(variable, variable=="V5", "Control6")) %>% 
  mutate(variable = replace(variable, variable=="V6", "Control10")) %>% 
  mutate(variable = replace(variable, variable=="V8", "Control6")) %>% 
  mutate(variable = replace(variable, variable=="V9", "Control10")) %>% 
  mutate(variable = replace(variable, variable=="V11", "Control6")) %>% 
  mutate(variable = replace(variable, variable=="V12", "Control10")) %>% 
  mutate(variable = replace(variable, variable=="V14", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V15", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V17", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V18", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V20", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V21", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V23", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V24", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V26", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V27", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V29", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V30", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V32", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V33", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V35", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V36", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V38", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V39", "Treatment10")) %>% 
  mutate(variable = replace(variable, variable=="V41", "Treatment6")) %>% 
  mutate(variable = replace(variable, variable=="V42", "Treatment10")) %>% 
  mutate(Treat = replace(Treat, variable=="Control6", "Control"),  Generation  = replace(Generation, variable=="Control6", "6")) %>% 
  mutate(Treat = replace(Treat, variable=="Control10", "Control"),  Generation  = replace(Generation, variable=="Control10", "10"))  %>% 
  mutate(Treat= replace(Treat, variable=="Treatment6", "Treatment"),  Generation  = replace(Generation, variable=="Treatment6", "6"))  %>% 
  mutate(Treat = replace(Treat, variable=="Treatment10", "Treatment"),  Generation  = replace(Generation, variable=="Treatment10", "10")) 

## Adding coverage of each
cov<-select(genoF_cov,2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33,35,36,38,39,41,42)
setDT(cov)
cov<-melt(cov)
MgenoF$Coverage <- cov$value

## Adding which beaker
beaker<-c(rep("C1", 2*nrow(genoF)), rep("C2", 2*nrow(genoF)),rep("C3", 2*nrow(genoF)),rep("C4", 2*nrow(genoF)),rep("T1", 2*nrow(genoF)),rep("T2", 2*nrow(genoF)),rep("T3", 2*nrow(genoF)),rep("T4", 2*nrow(genoF)),rep("T5", 2*nrow(genoF)),rep("T6", 2*nrow(genoF)),rep("T7", 2*nrow(genoF)),rep("T8", 2*nrow(genoF)),rep("T9", 2*nrow(genoF)),rep("T10", 2*nrow(genoF)))
MgenoF$Beaker<-beaker

## adjust data types and save
MgenoF$Generation <- as.numeric(MgenoF$Generation)
MgenoF$Coverage <- as.numeric(MgenoF$Coverage)
MgenoF$Beaker <- as.factor(MgenoF$Beaker)
MgenoF$Treat <- as.factor(MgenoF$Treat)
MgenoF <- as.data.table(MgenoF)

saveRDS(MgenoF,'full_data_lmm.raw.sim.1ksel_1knon.empNe.selstart.RDS')
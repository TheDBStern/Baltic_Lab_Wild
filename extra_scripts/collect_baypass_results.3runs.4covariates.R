library(gtools)
library(data.table)
library(matrixStats)
library(dplyr)

####################
## get C2 statistic
####################
contrast_files1 <- mixedsort(dir('.', pattern = 'run1_summary_contrast.out', full.names = TRUE))
contrast_files2 <- mixedsort(dir('.', pattern = 'run2_summary_contrast.out', full.names = TRUE))
contrast_files3 <- mixedsort(dir('.', pattern = 'run3_summary_contrast.out', full.names = TRUE))

contrast.tables1 = lapply(
  contrast_files1,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
contrast.tables2 = lapply(
  contrast_files2,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
contrast.tables3 = lapply(
  contrast_files3,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
run1 <- do.call(rbind,contrast.tables1)
run2 <- do.call(rbind,contrast.tables2)
run3 <- do.call(rbind,contrast.tables3)

con1_C2 <- rowMeans(as.matrix(data.frame(filter(run1,CONTRAST==1)$C2_std,filter(run2,CONTRAST==1)$C2_std,filter(run3,CONTRAST==1)$C2_std)))
con2_C2 <- rowMeans(as.matrix(data.frame(filter(run1,CONTRAST==2)$C2_std,filter(run2,CONTRAST==2)$C2_std,filter(run3,CONTRAST==2)$C2_std)))
con3_C2 <- rowMeans(as.matrix(data.frame(filter(run1,CONTRAST==3)$C2_std,filter(run2,CONTRAST==3)$C2_std,filter(run3,CONTRAST==3)$C2_std)))
con4_C2 <- rowMeans(as.matrix(data.frame(filter(run1,CONTRAST==4)$C2_std,filter(run2,CONTRAST==4)$C2_std,filter(run3,CONTRAST==4)$C2_std)))
con5_C2 <- rowMeans(as.matrix(data.frame(filter(run1,CONTRAST==5)$C2_std,filter(run2,CONTRAST==5)$C2_std,filter(run3,CONTRAST==5)$C2_std)))

###################
## get standardized XtX statistic
###################
xtx_files1 <- mixedsort(dir('.', pattern = 'run1_summary_pi_xtx.out', full.names = TRUE))
xtx_files2 <- mixedsort(dir('.', pattern = 'run2_summary_pi_xtx.out', full.names = TRUE))
xtx_files3 <- mixedsort(dir('.', pattern = 'run3_summary_pi_xtx.out', full.names = TRUE))

xtx.tables1 = lapply(
  xtx_files1,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
xtx.tables2 = lapply(
  xtx_files2,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
xtx.tables3 = lapply(
  xtx_files3,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
run1 <- do.call(rbind,xtx.tables1)
run2 <- do.call(rbind,xtx.tables2)
run3 <- do.call(rbind,xtx.tables3)

XtX <- rowMeans(as.matrix(data.frame(run1$XtXst,run2$XtXst,run3$XtXst)))

###################
## get beta correlation coefficients
###################
beta_files1 <- mixedsort(dir('.', pattern = 'run1_summary_betai.out', full.names = TRUE))
beta_files2 <- mixedsort(dir('.', pattern = 'run2_summary_betai.out', full.names = TRUE))
beta_files3 <- mixedsort(dir('.', pattern = 'run3_summary_betai.out', full.names = TRUE))

beta.tables1 = lapply(
  beta_files1,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
beta.tables2 = lapply(
  beta_files2,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
beta.tables3 = lapply(
  beta_files3,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
run1 <- do.call(rbind,beta.tables1)
run2 <- do.call(rbind,beta.tables2)
run3 <- do.call(rbind,beta.tables3)

Beta <- rowMeans(as.matrix(data.frame(run1$M_Beta,run2$M_Beta,run3$M_Beta)))


snpdet <- read.table("../../../../wild.noPBE.snpdet",header=F)
#res <- data.frame(Transcript=snpdet[,1],Position=snpdet[,2], XtXSt=XtX,C2_All=con1_C2,C2_All_noKiel=con2_C2,C2_N_Baltic=con3_C2,C2_N_Baltic_noKiel=con4_C2,C2_N_Sea=con5_C2)  
res <- data.frame(Pseudo_Transcript=snpdet[,1],Pseudo_Position=snpdet[,2], Beta_surface_Baltic=Beta)  

#res <- res[order(res$CHR,res$BP),]
#res <- res[mixedorder(as.vector(res$CHR)),]
saveRDS(res,'../../std_beta.Baltic.mean.RDS')


library(gtools)
library(data.table)
library(matrixStats)

##
std_files1 <- mixedsort(dir('.', pattern = 'run1_summary_betai', full.names = TRUE))
std_files2 <- mixedsort(dir('.', pattern = 'run2_summary_betai', full.names = TRUE))
std_files3 <- mixedsort(dir('.', pattern = 'run3_summary_betai', full.names = TRUE))

std.tables1 = lapply(
  std_files1,
  function(x) 
  {
  	print(x)
  	table <- fread(x, header=T,data.table=F)
    table
  }
)
std.tables2 = lapply(
  std_files2,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
std.tables3 = lapply(
  std_files3,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
std.snp.res1 <- do.call(rbind,std.tables1)
std.snp.res2 <- do.call(rbind,std.tables2)
std.snp.res3 <- do.call(rbind,std.tables3)

std_BF_all <- as.matrix(data.frame(std.snp.res1[,7],std.snp.res2[,7],std.snp.res3[,7]))
std_BF <- rowMedians(std_BF_all)
rm(std.tables1,std.tables2,std.tables3,std.snp.res1,std.snp.res2,std.snp.res3)

##
aux_files1 <- mixedsort(dir('aux_run1', pattern = 'betai', full.names = TRUE))
aux_files2 <- mixedsort(dir('aux_run2', pattern = 'betai', full.names = TRUE))
aux_files3 <- mixedsort(dir('aux_run3', pattern = 'betai', full.names = TRUE))

aux.tables1 = lapply(
  aux_files1,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
aux.tables2 = lapply(
  aux_files2,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
aux.tables3 = lapply(
  aux_files3,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
aux.snp.res1 <- do.call(rbind,aux.tables1)
aux.snp.res2 <- do.call(rbind,aux.tables2)
aux.snp.res3 <- do.call(rbind,aux.tables3)

aux_BF_all <- as.matrix(data.frame(aux.snp.res1[,6],aux.snp.res2[,6],aux.snp.res3[,6]))
aux_BF <- rowMedians(aux_BF_all)
##
core_files1 <- mixedsort(dir('.', pattern = 'run1', full.names = TRUE))
core_files2 <- mixedsort(dir('.', pattern = 'run2', full.names = TRUE))
core_files3 <- mixedsort(dir('.', pattern = 'run3', full.names = TRUE))

core.tables1 = lapply(
  core_files1,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
core.tables2 = lapply(
  core_files2,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
core.tables3 = lapply(
  core_files3,
  function(x) 
  {
    print(x)
    table <- fread(x, header=T,data.table=F)
    table
  }
)
core.snp.res1 <- do.call(rbind,core.tables1)
core.snp.res2 <- do.call(rbind,core.tables2)
core.snp.res3 <- do.call(rbind,core.tables3)

XtX_all <- as.matrix(data.frame(core.snp.res1[,6],core.snp.res2[,6],core.snp.res3[,6]))
XtX <- rowMeans(XtX_all)


snpdet <- read.table("../../../wild.noPBE.snpdet",header=F)
res <- data.frame(Transcript=snpdet[,1],Position=snpdet[,2], std_eBP=std_BF)  
#res <- res[order(res$CHR,res$BP),]
#res <- res[mixedorder(as.vector(res$CHR)),]
saveRDS(res,'../../../std_eBP.6pops.medians.RDS')


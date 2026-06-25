##
### ---------------
###
### ---------------

rm(list = ls())
options(stringsAsFactors = F)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require(bsseq))install.packages('bsseq')
install.packages('bsseq', repos= "https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages")
BiocManager::install("bsseq")
BiocManager::install("DMRcate")
library(bsseq)

if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("DMRcate")
BiocManager::install("bacon")


require(GEOquery)
require(Biobase)
library(BiocManager)
library("impute")
library(minfi)
library(minfiData)
library(sva)
library(GEOquery)
library(ChAMP)
library(bacon)
library(knitr)
library(limma)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(methylationArrayAnalysis)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)  ##470M
head(ann450k)
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) ##1GB
head(ann850k)

targetsall <- read.metharray.sheet("GSE111629/idat", pattern="567allsamplesheet.csv") 
targetsPD <- read.metharray.sheet("GSE111629/idat", pattern="331PDsamplesheet.csv")
targetsHC <- read.metharray.sheet("GSE111629/idat", pattern="236HCsamplesheet.csv")
targetsall
targetsPD
targetsH
rm (targetsall)
rm (targetsPD)
rm (targetsHC)
targetsall <- read.metharray.sheet("GSE111629/idat", pattern="566allsamplesheet.csv") 


rgSetall <- read.metharray.exp("GSE111629/idat", targets=targetsall)
rgSetPD <- read.metharray.exp("GSE111629/idat", targets=targetsPD)
rgSetHC <- read.metharray.exp("GSE111629/idat", targets=targetsHC)
rgSetall <- read.metharray.exp("GSE111629/idat", targets=targetsall)
save(rgSetall,file = 'GSE111629_minfi_rgSetall566.Rdata')
load(file = 'GSE111629_minfi_rgSetall566.Rdata')



detP <- detectionP(rgSetall)
head(detP)
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targetsall$Sample_Name)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
failed <- colMeans(detP)>0.01 
rm (failed)
rm (detP)
rm (pal)


getSex(object = rgSetall, cutoff = -2)
class(rgSetall)
Msetall <- preprocessRaw(rgSetall)
GMsetall <- mapToGenome(Msetall)

PreSex <- getSex(object = GMsetall, cutoff = -2)
class(PreSex) #DFrame

pd <- pData(rgSetall)
class(pd) #DFrame
install.packages(do)
library(do)
pdsex <- Replace(pd$Sex,"Female", "F") 
pdsex <- Replace(pdsex,"Male", "M")
pdsex == PreSex$predictedSex


library(FlowSorted.Blood.450k)
cellCounts <- estimateCellCounts(rgSetall)
save(cellCounts,file = 'GSE111629_minfi_cellCounts.csv')



myLoadHC <-champ.load("GSE111629/idat",
           method="ChAMP",
           methValue="B",
           autoimpute=TRUE,
           filterDetP=TRUE,
           ProbeCutoff=0,
           SampleCutoff=0.1,
           detPcut=0.01,
           filterBeads=TRUE,
           beadCutoff=0.05,
           filterNoCG=TRUE,
           filterSNPs=TRUE,
           population=NULL,
           filterMultiHit=TRUE,
           filterXY=TRUE,
           force=FALSE,
           arraytype="450K")

save(myLoadHC,file = 'GSE111629_champ_myLoadHC236.Rdata')
load(file = 'GSE111629_champ_myLoadHC236.Rdata')

champ.QC(beta = myLoadHC$beta,
         pheno=myLoadHC$pd$Sex,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir="./CHAMP_QCimages_236HC/")

myNormHC <- champ.norm(beta=myLoadHC$beta,
                       rgSet=myLoadHC$rgSet,
                       mset=myLoadHC$mset,
                       resultsDir="./CHAMP_Normalization_236HC/",
                       method="BMIQ",
                       plotBMIQ=FALSE,
                       arraytype="450K",
                       cores=3)
class(myNormHC)   #"matrix" "array" 
dim(myNormHC)  #[1] 407788    236
pD=myLoadHC$pd
pD$Age<-as.numeric(as.character(pD$Age))
pD$Slide<-as.character(as.numeric(pD$Slide))
save(myNormHC,pD,file = 'GSE111629_champ_myNormHC236.Rdata')
load(file = 'GSE111629_champ_myNormHC236.Rdata')


champ.SVD(beta = myNormPD,
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_330PD/")

class(myNormPD) %>% length()   #equals 2
class(myNormPD) == "data.frame"   #So it will cause problems:  the condition has length > 1 and only the first element will be used
## u can code this champ.SVD(beta=myNorm %>% as.data.frame(), pd=myLoad$pd) to avoid it

champ.SVD(beta = myNormHC %>% as.data.frame(),
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_236HC/")


?champ.runCombat
myCombatHC <- champ.runCombat(beta=myNormHC,
                              pd=pD,
                              variablename="Sex",
                              batchname=c("Slide"),
                              logitTrans=TRUE)

champ.SVD(beta = myCombatHC %>% as.data.frame(),
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_236HC_AftercombatSlide/")

myCombatHC <- champ.runCombat(beta=myCombatHC,
                              pd=pD,
                              variablename="Sex",
                              batchname=c("Array"),
                              logitTrans=TRUE)

champ.SVD(beta = myCombatHC %>% as.data.frame(),
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_236HC_AftercombatSlideArray/")

save(myCombatHC,pD,file = 'GSE111629_champ_myCombatHC236.Rdata')
load(file = 'GSE111629_champ_myCombatHC236.Rdata')


?champ.refbase
myRefbaseHC <- champ.refbase(beta = myCombatHC,
              arraytype="450K")
class(myRefbaseHC)  #[1] "list"
class(myRefbaseHC$CorrectedBeta)   #[1] "matrix" "array" 
class(myRefbaseHC$CellFraction)  #[1] "matrix" "array" 
save(myRefbaseHC,pD,file = 'GSE111629_champ_myRefbaseHC236.Rdata')
load(file = 'GSE111629_champ_myRefbaseHC236.Rdata')


?champ.DMP
myDMPHC0.05 <- champ.DMP(beta = myRefbaseHC$CorrectedBeta,
                     pheno = pD$Sex,
                     compare.group = NULL,
                     adjPVal = 0.05,
                     adjust.method = "BH",
                     arraytype = "450K")
head(myDMPHC0.05[[1]]) 
dim(myDMPHC0.05$F_to_M)
save(myDMPHC0.05,file = 'GSE111629_champ_myDMPHC0.05.Rdata')
load(file = 'GSE111629_champ_myDMPHC0.05.Rdata')
write.table(myDMPHC0.05, file = "myDMPHC0.05.csv", sep = ",")

myDMP <- champ.DMP(beta = myNormHC,
                   pheno = pD$Sex,
                   compare.group = NULL,
                   adjPVal = 0.05,
                   adjust.method = "BH",
                   arraytype = "450K")
head(myDMP[[1]])
dim(myDMP$F_to_M)
myDMP <- champ.DMP(beta = myCombatHC,
                   pheno = pD$Sex,
                   compare.group = NULL,
                   adjPVal = 0.05,
                   adjust.method = "BH",
                   arraytype = "450K")
head(myDMP[[1]])
dim(myDMP$F_to_M)
rm(myDMP)


myDMPHC1 <- champ.DMP(beta = myRefbaseHC$CorrectedBeta,
                         pheno = pD$Sex,
                         compare.group = NULL,
                         adjPVal = 1,
                         adjust.method = "BH",
                         arraytype = "450K")
head(myDMPHC1[[1]]) 
dim(myDMPHC1$F_to_M)
save(myDMPHC1,file = 'GSE111629_champ_myDMPHC1.Rdata')
load(file = 'GSE111629_champ_myDMPHC1.Rdata')
write.table(myDMPHC1, file = "myDMPHC1.csv", sep = ",")

library(qqman)
qq(myDMPHC1$F_to_M$P.Value)
p_value=myDMPHC1$F_to_M$P.Value
z = qnorm(p_value/ 2)
lambda = median(z^2, na.rm = TRUE) / 0.454
lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)
lambda = round(median(z^2, na.rm = TRUE) / qchisq(0.5, 1), 3)
lambda = round(median(z^2, na.rm = TRUE) / 0.456, 3)

qq(myDMPHC1$F_to_M$adj.P.Val)  
p_value=myDMPHC1$F_to_M$adj.P.Val
z = qnorm(p_value/ 2)
lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)
lambda = round(median(z^2, na.rm = TRUE) / qchisq(0.5, 1), 3)
lambda = round(median(z^2, na.rm = TRUE) / 0.456, 3)

Estimate = logFC
StdErr = Estimate / tvalue = logFC / tvalue
Estimate = logFC
StdErr = logFC / tvalue
Estimate = myDMPHC1$F_to_M$logFC
StdErr = Estimate/myDMPHC1$F_to_M$t
## bacon analysis
bc <- bacon(
  teststatistics = NULL,
  effectsizes =  Estimate,
  standarderrors = StdErr,
  na.exclude = TRUE
)

## inflation factor
print("lambda.bacon")
print(inflation(bc))  #1.231635 

## data$tvalue <- data$Estimate / data$StdErr
## data$zvalue <- qnorm(pt(data$tvalue, df))
## data$chisq <- (data$zvalue) ^ 2
Estimate = myDMPHC1$F_to_M$logFC
StdErr = Estimate/myDMPHC1$F_to_M$t
zvalue = qnorm(myDMPHC1$F_to_M$P.Value / 2)
inflationFactor = median(zvalue^2,na.rm = TRUE) / qchisq(0.5, 1)
print("lambda")
print(inflationFactor)

# genome-wide sig cpgs
sig <- ifelse(myDMPHC1$F_to_M$P.Value < 2.4e-7, 1, 0)
table(sig)

### 2. bacon analysis
bc <- bacon(
  teststatistics = NULL,
  effectsizes =  Estimate,
  standarderrors = StdErr,
  na.exclude = TRUE
)

# inflation factor
print("lambda.bacon")
print(inflation(bc))

### 3. Create final dataset
myDMPHC_bacon <- data.frame(
  myDMPHC1,
  Estimate.bacon = bacon::es(bc),
  StdErr.bacon = bacon::se(bc),
  pValue.bacon = pval(bc),
  fdr.bacon = p.adjust(pval(bc), method = "fdr"),
  stringsAsFactors = FALSE)

class(myDMPHC_bacon)  # "data.frame"
keep <- myDMPHC_bacon$F_to_M.adj.P.Val < 0.05  
myDMPHC0.05_bacon = myDMPHC_bacon[keep,]
write.table(myDMPHC0.05_bacon, file = "myDMPHC0.05_bacon.csv", sep = ",")








?champ.DMR
myDMRHC0.05 <- champ.DMR(beta=myRefbaseHC$CorrectedBeta,pheno=pD$Sex,method="Bumphunter")
write.table(myDMRHC0.05, file = "myDMRHC0.05.csv", sep = ",")





myLoad <-champ.load("GSE111629/idat",
                      method="ChAMP",
                      methValue="B",
                      autoimpute=TRUE,
                      filterDetP=TRUE,
                      ProbeCutoff=0,
                      SampleCutoff=0.1,
                      detPcut=0.01,
                      filterBeads=TRUE,
                      beadCutoff=0.05,
                      filterNoCG=TRUE,
                      filterSNPs=TRUE,
                      population=NULL,
                      filterMultiHit=TRUE,
                      filterXY=TRUE,
                      force=FALSE,
                      arraytype="450K")

save(myLoad,file = 'GSE111629_myLoad567.Rdata')
load(file = 'GSE111629_myLoad567.Rdata')


myLoadPD <-champ.load("GSE111629/idat",
                    method="ChAMP",
                    methValue="B",
                    autoimpute=TRUE,
                    filterDetP=TRUE,
                    ProbeCutoff=0,
                    SampleCutoff=0.1,
                    detPcut=0.01,
                    filterBeads=TRUE,
                    beadCutoff=0.05,
                    filterNoCG=TRUE,
                    filterSNPs=TRUE,
                    population=NULL,
                    filterMultiHit=TRUE,
                    filterXY=TRUE,
                    force=FALSE,
                    arraytype="450K")








untar("GSE111629_RAW.tar", exdir = "GSE111629/idat")
head(list.files("GSE111629/idat", pattern = "idat"))
idatFiles <- list.files("GSE111629/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

rgSet <- read.metharray.exp("GSE111629/idat")
rgSet
pData(rgSet)  ##DataFrame with 40 rows and 0 columns
head(sampleNames(rgSet))
class(rgSet)  ##[1] "RGChannelSet", attr(,"package"), [1] "minfi"
save(rgSet,file = 'GSE111629_minfi_rgSet.Rdata')
load(file = 'GSE111629_minfi_rgSet.Rdata')

dataDirectory <- file("path/to/PDSex/blood_111629_571/ALL571")
list.files(dataDirectory, recursive = TRUE)
rm(dataDirectory)
targets <- read.metharray.sheet(dataDirectory, pattern="571sampleSheet.csv")
targets <- read.metharray.sheet(pattern="571samplesheet.csv")
targets <- read.csv(file = "571samplesheet.csv")

sampleNames(rgSet) <- targets$Sample_Name
rgSet

##



detP <- detectionP(rgSet)
head(detP)
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
failed <- colMeans(detP)>0.01
rm (failed)
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")
qcReport(rgSet, sampNames = targets$ID,sampGroups = targets$Sample_Group,
         pdf="qcReport.pdf")

keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet

targets <- targets[keep,]
targets[,1:5]

detP <- detP[,keep]
dim(detP)

mSetSq <- preprocessQuantile(rgSet) 
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))


par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",  
        col=pal[factor(targets$Sample_Source)])
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       bg="white", cex=0.7)

par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal, 
       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(3,4))
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")


detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt


keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep) #TRUE：472927；FALSE：11608。472927+11608=484535
mSetSqFlt <- mSetSqFlt[keep,]

mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE) 
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)  #TRUE：429397；FALSE：26527。429397+26527=455924
mSetSqFlt <- mSetSqFlt[keep,] 
mSetSqFlt

par(mfrow=c(1,2))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], cex=0.8)
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.45, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Source)])
legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")
par(mfrow=c(1,3))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Source)], dim=c(1,3))
legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Source)], dim=c(2,3))
legend("topright", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Source)], dim=c(3,4))
legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")


bVals <- getBeta(mSetSqFlt)
class(mSetSqFlt) #"GenomicRatioSet" attr(,"package") [1] "minfi"
dim(bVals) #[1] 429397     10，
class(bVals)
library(ChAMP)
par(mfrow=c(1,1))
champ.SVD(beta = bVals , 
          rgSet=NULL,
          pd=targets,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages/")


mVals <- getM(mSetSqFlt)
head(mVals[,1:5])
bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))










info=read.table("group.txt",sep="\t",header=T)
library(data.table)
b=info
rownames(b)=b[,1]
a=fread("data.txt",data.table = F )
a[1:4,1:4]
rownames(a)=a[,1]
a=a[,-1]
beta=as.matrix(a)
beta=impute.knn(beta)
betaData=beta$data
betaData=betaData+0.00001
a=betaData
a[1:4,1:4]
identical(colnames(a),rownames(b))

library(ChAMP)
myLoad=champ.filter(beta = a,pd = b)
myLoad
save(myLoad,file = 'step1-output.Rdata')


if(F){
  require(GEOquery)
  require(Biobase)
  eset <- getGEO("GSE144858",destdir = './',AnnotGPL = T,getGPL = F)
  beta.m <- exprs(eset[[1]])
  pD.all <- pData(eset[[1]])
  pd=subset(pD.all, description=="Alzheimer's disease sample")
  pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
  pD <- pD.all[, c(36, 40)]
  head(pD) 
  names(pD)[c(1,2)] <- c("age", "sex")
  pD$group <- sub("^diagnosis: ", "", pD$group)
  pD$sex <- sub("^Sex: ", "", pD$sex)
  pD1 <- t(pD)
  beta.m=subset(beta.m0, select=c(GSM4299203,GSM4299204,GSM4299208,GSM4299211,GSM4299212,GSM4299220,GSM4299221,GSM4299222,GSM4299223,GSM4299225,GSM4299227,GSM4299228,GSM4299230,GSM4299234,GSM4299235,GSM4299239,GSM4299243,GSM4299349,GSM4299350,GSM4299354,GSM4299358,GSM4299359,GSM4299360,GSM4299365,GSM4299366,GSM4299368,GSM4299369,GSM4299371,GSM4299373,GSM4299379,GSM4299380,GSM4299382,GSM4299383,GSM4299388,GSM4299390,GSM4299391,GSM4299392,GSM4299393,GSM4299400,GSM4299403,GSM4299408,GSM4299410,GSM4299411,GSM4299417,GSM4299431,GSM4299432,GSM4299433,GSM4299436,GSM4299438,GSM4299453,GSM4299455,GSM4299456,GSM4299457,GSM4299470,GSM4299472,GSM4299473,GSM4299476,GSM4299477,GSM4299485,GSM4299487,GSM4299490,GSM4299495,GSM4299551,GSM4299555,GSM4299556,GSM4299559,GSM4299569,GSM4299571,GSM4299572,GSM4299574,GSM4299575,GSM4299577,GSM4299578,GSM4299579,GSM4299582,GSM4299586,GSM4299588,GSM4299620,GSM4299622,GSM4299627,GSM4299636,GSM4299638,GSM4299639,GSM4299640,GSM4299641,GSM4299644,GSM4299650,GSM4299653,GSM4299655,GSM4299657,GSM4299662,GSM4299663,GSM4299664))

  library(ChAMP)
  beta=beta.m
  beta=impute.knn(beta)
  betaData=beta$data
  betaData=betaData+0.00001
  beta.m=betaData
  beta.m[1:4,1:4]
  identical(colnames(beta.m),rownames(pD))
  
  myLoad=champ.filter(beta = beta.m ,pd = pD)
  myLoad
  save(myLoad,file = 'step1-output.Rdata')
}


champ.QC(beta = myLoad$beta,
         pheno=myLoad$pd$sex,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir="./CHAMP_QCimages/")



if(F){
  myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=3)
  dim(myNorm) 
  pD=myLoad$pd
  save(myNorm,pD,file = 'step2-champ_myNorm.Rdata')
}
load(file = 'step2-champ_myNorm.Rdata')
beta.m=myNorm

champ.SVD(beta=myNorm,pd=myLoad$pd)
pD$age<-as.numeric(as.character(pD$age))
myLoad$pd$age<-as.numeric(as.character(myLoad$pd$age))

##DMPs
myDMP <- champ.DMP(beta = myNorm,
          pheno = myLoad$pd$sex,
          compare.group = NULL,
          adjPVal = 0.05,
          adjust.method = "BH",
          arraytype = "450K")
head(myDMP[[1]])
save(myDMP,file = 'step3-output-myDMP.Rdata')

myDMP0.1 <- champ.DMP(beta = myNorm,
                   pheno = myLoad$pd$sex,
                   compare.group = NULL,
                   adjPVal = 0.1,
                   adjust.method = "BH",
                   arraytype = "450K")

DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=myLoad$pd$sex)

DMP.GUI(DMP=myDMP[[1]],
       beta=myNorm,
       pheno=myLoad$pd$sex,
       cutgroupnumber=4)
write.table(myDMP0.1, file = "myDMP0.1-AD93.csv", sep = ",")
write.table(myDMP, file = "myDMP369-AD93.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

myDMR2 <- champ.DMR(beta=myNorm,
          pheno=myLoad$pd$sex,
          compare.group=NULL,
          arraytype="450K",
          method = "Bumphunter",
          minProbes=7,
          adjPvalDmr=0.05,
          cores=3,
          ## following parameters are specifically for Bumphunter method.
          maxGap=300,
          cutoff=NULL,
          pickCutoff=TRUE,
          smooth=TRUE,
          smoothFunction=loessByCluster,
          useWeights=FALSE,
          permutations=NULL,
          B=250,
          nullMethod="bootstrap")
save(myDMR,file = 'step4-output-myDMR.Rdata')
write.table(myDMR, file = "myDMR-AD93.csv", sep = ",")

write.table(myDMR1, file = "myDMR1-AD93.csv", sep = ",")
write.table(myDMR2, file = "myDMR2-AD93.csv", sep = ",")

DMR.GUI(DMR=myDMR,
        beta=myNorm,
        pheno=myLoad$pd$sex,
        runDMP=TRUE,
        compare.group=NULL,
        arraytype="450K")


myGSEA <- champ.GSEA(beta=myNorm,
           DMP=myDMP[[1]],
           DMR=myDMR,
           CpGlist=NULL,
           Genelist=NULL,
           pheno=myLoad$pd$sex,
           method="fisher",
           arraytype="450K",
           Rplot=TRUE,
           adjPval=0.05,
           cores=1)
write.table(myGSEA$DMP, file = "myGSEA_DMP-AD93.csv", sep = ",")
write.table(myGSEA$DMR, file = "myGSEA_DMR-AD93.csv", sep = ",")
myGSEAfisher <- myGSEA
rm(myGSEA)

myGSEAgometh <- champ.GSEA(beta=myNorm,
                     DMP=myDMP[[1]],
                     DMR=myDMR,
                     CpGlist=NULL,
                     Genelist=NULL,
                     pheno=myLoad$pd$sex,
                     method="gometh",
                     arraytype="450K",
                     Rplot=TRUE,
                     adjPval=0.05,
                     cores=1)


myGSEAeb <- champ.ebGSEA(beta=myNorm, pheno=myLoad$pd$sex, minN=5, adjPval=0.05, arraytype="450K", cores=1)


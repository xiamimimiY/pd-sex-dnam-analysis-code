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

targetsall <- read.metharray.sheet("GSE145361/idat", pattern="1889allsamplesheet.csv") 
rm (targetsall)
class(targetsall)
del <- which(targetsall$Basename=="character(0)")
targetsall <- targetsall[-del,]  ##1889-298=1591

rgSetall <- read.metharray.exp("GSE145361/idat",targets=targetsall)
rgSet <- read.metharray.exp("GSE145361/idat")
dim(rgSet) #[1] 622399   1590
rgSetall = rgSet
rm (rgSet)
class(rgSetall)  #[1] "RGChannelSet" attr(,"package")  [1] "minfi"

rgSetPD <- read.metharray.exp("GSE111629/idat", targets=targetsPD)
rgSetHC <- read.metharray.exp("GSE111629/idat", targets=targetsHC)
rgSetall <- read.metharray.exp("GSE111629/idat", targets=targetsall)
save(rgSetall,file = 'GSE111629_minfi_rgSetall566.Rdata')
load(file = 'GSE111629_minfi_rgSetall566.Rdata')



detP <- detectionP(rgSetall)
dim(detP)  #[1] 485512   1590
head(detP)
?ColorBrewer
install.packages('RColorBrewer')
library(RColorBrewer)
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targetsall$Sample_Name)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
failed <- colMeans(detP)>0.01
keep <- colMeans(detP) < 0.01
rm (failed)
rm (detP)
rm (pal)
rm (keep)

getSex(object = rgSetall, cutoff = -2)
class(rgSetall)
Msetall <- preprocessRaw(rgSetall)
dim(Msetall)  #[1] 485512   1590
GMsetall <- mapToGenome(Msetall)
dim(GMsetall)  #[1] 485512   1590

PreSex <- getSex(object = GMsetall, cutoff = -2)
class(PreSex) #DFrame

targetsall$Sex == PreSex$predictedSex
all(targetsall$Sex == PreSex$predictedSex)
PreSex$predictedSex == X1590samplesheet$Sex
all(PreSex$predictedSex == X1590samplesheet$Sex)
rm (X1590samplesheet)


library(FlowSorted.Blood.450k)
cellCounts <- estimateCellCounts(rgSetall)
save(cellCounts,file = 'GSE111629_minfi_cellCounts.csv')



myLoadPD <-champ.load("GSE145361/idat",
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


champ.QC(beta = myLoadPD$beta,
         pheno=myLoadPD$pd$Sex,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir="./CHAMP_QCimages_813PD/")

myNormPD <- champ.norm(beta=myLoadPD$beta,
                       rgSet=myLoadPD$rgSet,
                       mset=myLoadPD$mset,
                       resultsDir="./CHAMP_Normalization_813PD/",
                       method="BMIQ",
                       plotBMIQ=FALSE,
                       arraytype="450K",
                       cores=3)
class(myNormPD)   #"matrix" "array" 
dim(myNormPD)  #[1] 404987    813
pD=myLoadPD$pd
pD$Slide<-as.character(as.numeric(pD$Slide))
save(myNormPD,pD,file = 'GSE145361_champ_myNormPD813.Rdata')   #12：59-
load(file = 'GSE145361_champ_myNormPD813.Rdata')


champ.SVD(beta = myNormPD,
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_813PD/")

class(myNormPD) %>% length()   #equals 2
class(myNormPD) == "data.frame"   #So it will cause problems:  the condition has length > 1 and only the first element will be used
## u can code this champ.SVD(beta=myNorm %>% as.data.frame(), pd=myLoad$pd) to avoid it

champ.SVD(beta = myNormPD %>% as.data.frame(),
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_813PD/")


?champ.runCombat
myCombatPD <- champ.runCombat(beta=myNormPD,
                              pd=pD,
                              variablename="Sex",
                              batchname=c("Slide"),
                              logitTrans=TRUE)

champ.SVD(beta = myCombatPD %>% as.data.frame(),
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_813PD_AftercombatSlide/")

myCombatPD <- champ.runCombat(beta=myCombatPD,
                              pd=pD,
                              variablename="Sex",
                              batchname=c("Array"),
                              logitTrans=TRUE)
save(myCombatPD,pD,file = 'GSE145361_champ_myCombatPD813.Rdata')
load(file = 'GSE145361_champ_myCombatPD813.Rdata')


champ.SVD(beta = myCombatPD %>% as.data.frame(),
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_813PD_AftercombatSlideArray/")


?champ.refbase
myRefbasePD <- champ.refbase(beta = myCombatPD,
              arraytype="450K")
class(myRefbasePD)  #[1] "list"
class(myRefbasePD$CorrectedBeta)   #[1] "matrix" "array" 
class(myRefbasePD$CellFraction)  #[1] "matrix" "array" 
class(myNormPD)    #[1] "matrix" "array"
save(myRefbasePD,pD,file = 'GSE145361_champ_myRefbasePD813.Rdata')   
load(file = 'GSE145361_champ_myRefbasePD813.Rdata')


?champ.DMP
myDMPPD0.05_813 <- champ.DMP(beta = myRefbasePD$CorrectedBeta,
                     pheno = pD$Sex,
                     compare.group = NULL,
                     adjPVal = 0.05,
                     adjust.method = "BH",
                     arraytype = "450K")
head(myDMPPD0.05_813[[1]]) 
dim(myDMPPD0.05_813$M_to_F)
save(myDMPPD0.05_813,file = 'GSE145361_champ_myDMPPD0.05_813.Rdata')
load(file = 'GSE145361_champ_myDMPPD0.05_813.Rdata')
write.table(myDMPPD0.05_813, file = "myDMPPD0.05_813.csv", sep = ",")

myDMP <- champ.DMP(beta = myNormPD,
                   pheno = pD$Sex,
                   compare.group = NULL,
                   adjPVal = 0.05,
                   adjust.method = "BH",
                   arraytype = "450K")
head(myDMP[[1]])
dim(myDMP$M_to_F)
myDMP <- champ.DMP(beta = myCombatPD,
                   pheno = pD$Sex,
                   compare.group = NULL,
                   adjPVal = 0.05,
                   adjust.method = "BH",
                   arraytype = "450K")
head(myDMP[[1]])
dim(myDMP$M_to_F)
rm(myDMP)


myDMPPD1_813 <- champ.DMP(beta = myRefbasePD$CorrectedBeta,
                         pheno = pD$Sex,
                         compare.group = NULL,
                         adjPVal = 1,
                         adjust.method = "BH",
                         arraytype = "450K")
head(myDMPPD1_813[[1]]) 
dim(myDMPPD1_813$M_to_F)
save(myDMPPD1_813,file = 'GSE145361_champ_myDMPPD1_813.Rdata')
load(file = 'GSE145361_champ_myDMPPD1_813.Rdata')


library(qqman)
qq(myDMPPD1_813$M_to_F$P.Value)
p_value=myDMPPD1_813$M_to_F$P.Value
z = qnorm(p_value/ 2)
lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)
lambda = round(median(z^2, na.rm = TRUE) / qchisq(0.5, 1), 3)
lambda = round(median(z^2, na.rm = TRUE) / 0.456, 3)



## data$tvalue <- data$Estimate / data$StdErr
## data$zvalue <- qnorm(pt(data$tvalue, df))
## data$chisq <- (data$zvalue) ^ 2
Estimate = myDMPPD1_813$M_to_F$logFC
StdErr = Estimate/myDMPPD1_813$M_to_F$t
zvalue = qnorm(myDMPPD1_813$M_to_F$P.Value / 2)
inflationFactor = median(zvalue^2,na.rm = TRUE) / qchisq(0.5, 1)
print("lambda")
print(inflationFactor)

# genome-wide sig cpgs
sig <- ifelse(myDMPPD1_813$M_to_F$P.Value < 2.4e-7, 1, 0)
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
myDMPPD1_813_bacon <- data.frame(
  myDMPPD1_813,
  Estimate.bacon = bacon::es(bc),
  StdErr.bacon = bacon::se(bc),
  pValue.bacon = pval(bc),
  fdr.bacon = p.adjust(pval(bc), method = "fdr"),
  stringsAsFactors = FALSE)

class(myDMPPD1_813_bacon)  # "data.frame"
keep <- myDMPPD1_813_bacon$M_to_F.adj.P.Val < 0.05  
myDMPPD813_0.05_bacon = myDMPPD1_813_bacon[keep,]
write.table(myDMPPD813_0.05_bacon, file = "myDMPPD813_0.05_bacon.csv", sep = ",")




?champ.DMR
myDMRPD0.05_813 <- champ.DMR(beta=myRefbasePD$CorrectedBeta,pheno=pD$Sex,method="Bumphunter")
write.table(myDMRPD0.05_813, file = "myDMRPD0.05_813.csv", sep = ",")


?DMR.GUI
DMR.GUI(DMR=myDMRPD0.05_813,
       beta=myRefbasePD$CorrectedBeta,
       pheno=pD$Sex,
       runDMP=TRUE,
       compare.group=NULL,
       arraytype="450K")













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


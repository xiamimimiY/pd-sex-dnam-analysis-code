##
### ---------------
###
###
### ---------------

rm(list = ls())
options(stringsAsFactors = F)

require(GEOquery)
require(Biobase)
library("impute")
library(minfi)
library(GEOquery)
library(ChAMP)

untar("GSE111629_RAW.tar", exdir = "GSE111629/idat")
head(list.files("GSE111629/idat", pattern = "idat"))

## Find CSV Success
## Reading CSV File
## Your pd file contains NO Array(Sentrix_Position) information.
## Your pd file contains NO Slide(Sentrix_ID) information.
## There is NO Pool_ID in your pd file.
## There is NO Sample_Plate in your pd file.
## There is NO Sample_Well in your pd file.
## Error in champ.import(directory, arraytype = arraytype) : 
## Error Match between pd file and Green Channel IDAT file.


myload <- champ.load(directory = getwd(),
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




idatFiles <- list.files("GSE111629/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

library("minfi")
rgSet <- read.metharray.exp("GSE111629/idat")
rgSet
pData(rgSet)  ##DataFrame with 40 rows and 0 columns
head(sampleNames(rgSet))
class(rgSet)  ##[1] "RGChannelSet", attr(,"package"), [1] "minfi"
save(rgSet,file = 'GSE111629_minfi_rgSet.Rdata')
load(file = 'GSE111629_minfi_rgSet.Rdata')

geoMat <- getGEO("GSE111629")
pD.all <- pData(geoMat[[1]])
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
head(pD)  
names(pD)[c(3,4)] <- c("group", "sex")
pD$group <- sub("^diagnosis: ", "", pD$group)
pD$sex <- sub("^Sex: ", "", pD$sex)

sampleNames(rgSet) <- sub(".*_5", "5", sampleNames(rgSet))
rownames(pD) <- pD$title
pD <- pD[sampleNames(rgSet),]
pData(rgSet) <- pD
rgSet




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


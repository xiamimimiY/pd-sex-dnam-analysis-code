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

##——————————————————2023-3-22————————————————————————————————————————————————————     
untar("GSE111223_RAW.tar", exdir = "GSE111223/idat")
head(list.files("GSE111223/idat", pattern = "idat"))     
     
targetsall <- read.metharray.sheet("GSE111223/idat", pattern="259allsamplesheet.csv") 

rgSetall <- read.metharray.exp("GSE111223/idat",targets=targetsall)
dim(rgSetall) #[1] 622399    259
class(rgSetall)  #[1] "RGChannelSet" attr(,"package")  [1] "minfi"

rgSetPD <- read.metharray.exp("GSE111629/idat", targets=targetsPD)
rgSetHC <- read.metharray.exp("GSE111629/idat", targets=targetsHC)
rgSetall <- read.metharray.exp("GSE111629/idat", targets=targetsall)
save(rgSetall,file = 'GSE111629_minfi_rgSetall566.Rdata')
load(file = 'GSE111629_minfi_rgSetall566.Rdata')



detP <- detectionP(rgSetall)
dim(detP)  #[1] 485512    259
head(detP)
?ColorBrewer
library(RColorBrewer)
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targetsall$Sample_Name)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
failed <- colMeans(detP)>0.01
keep <- colMeans(detP) < 0.01
print(max(colMeans(detP)))
print(min(colMeans(detP)))
rm (failed)
rm (detP)
rm (pal)
rm (keep)

getSex(object = rgSetall, cutoff = -2)
class(rgSetall)
Msetall <- preprocessRaw(rgSetall)
dim(Msetall)  #[1] 485512   259
GMsetall <- mapToGenome(Msetall)
dim(GMsetall)  #[1] 485512  259

PreSex <- getSex(object = GMsetall, cutoff = -2)
class(PreSex) #DFrame

targetsall$sex == PreSex$predictedSex
all(targetsall$Sex == PreSex$predictedSex)


library(FlowSorted.Blood.450k)
cellCounts <- estimateCellCounts(rgSetall)
save(cellCounts,file = 'GSE111629_minfi_cellCounts.csv')


library(data.table)
library(stringr)
untar("GSE111223_RAW.tar", exdir = "GSE111223/idat")
dir_path <- "GSE111223/idat"

idatFiles <- list.files(dir_path, pattern = ".idat.gz$", full.names = TRUE)

library(utils)
install.packages("utils")

for (i in seq_along(idatFiles)) {
  gunzip(idatFiles[i])
}


myLoadPD <-champ.load("GSE111223/idat",
           method="minfi",
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
         pheno=myLoadPD$pd$sex,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir="./CHAMP_QCimages_128PD/")

myNormPD <- champ.norm(beta=myLoadPD$beta,
                       rgSet=myLoadPD$rgSet,
                       mset=myLoadPD$mset,
                       resultsDir="./CHAMP_Normalization_128PD/",
                       method="BMIQ",
                       plotBMIQ=FALSE,
                       arraytype="450K",
                       cores=3)
class(myNormPD)   #"matrix" "array" 
dim(myNormPD)  #[1] 397218    128
pD=myLoadPD$pd
pD$Slide<-as.character(as.numeric(pD$Slide))
pD$age<-as.numeric(as.character(pD$age))
save(myNormPD,pD,file = 'GSE111223_champ_myNormPD128.Rdata')  
load(file = 'GSE111223_champ_myNormPD128.Rdata')


champ.SVD(beta = myNormPD,
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_128PD/")

class(myNormPD) %>% length()   #equals 2
class(myNormPD) == "data.frame"   #So it will cause problems:  the condition has length > 1 and only the first element will be used
## u can code this champ.SVD(beta=myNorm %>% as.data.frame(), pd=myLoad$pd) to avoid it

champ.SVD(beta = myNormPD %>% as.data.frame(),
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_128PD/")


?champ.runCombat
myCombatPD <- champ.runCombat(beta=myNormPD,
                              pd=pD,
                              variablename="sex",
                              batchname=c("Slide"),
                              logitTrans=TRUE)

champ.SVD(beta = myCombatPD %>% as.data.frame(),
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_128PD_AftercombatSlide/")

myCombatPD <- champ.runCombat(beta=myCombatPD,
                              pd=pD,
                              variablename="sex",
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
          resultsDir="./CHAMP_SVDimages_128PD_AftercombatSlideArray/")


?champ.refbase
myRefbasePD <- champ.refbase(beta = myCombatPD,
              arraytype="450K")
class(myRefbasePD)  #[1] "list"
class(myRefbasePD$CorrectedBeta)   #[1] "matrix" "array" 
class(myRefbasePD$CellFraction)  #[1] "matrix" "array" 
class(myNormPD)    #[1] "matrix" "array"
save(myRefbasePD,pD,file = 'GSE111223_champ_myRefbasePD128.Rdata')   
load(file = 'GSE111223_champ_myRefbasePD128.Rdata')


?champ.DMP
myDMPPD0.05_128 <- champ.DMP(beta = myCombatPD,
                     pheno = pD$sex,
                     compare.group = NULL,
                     adjPVal = 0.05,
                     adjust.method = "BH",
                     arraytype = "450K")
head(myDMPPD0.05_128[[1]]) 
dim(myDMPPD0.05_128$M_to_F)
save(myDMPPD0.05_128,file = 'GSE111223_champ_myDMPPD0.05_128.Rdata')
load(file = 'GSE111223_champ_myDMPPD0.05_128.Rdata')
write.table(myDMPPD0.05_128, file = "myDMPPD0.05_128.csv", sep = ",")

myDMP <- champ.DMP(beta = myNormPD,
                   pheno = pD$sex,
                   compare.group = NULL,
                   adjPVal = 0.05,
                   adjust.method = "BH",
                   arraytype = "450K")
head(myDMP[[1]])
dim(myDMP$M_to_F)
rm(myDMP)


myDMPPD1_128 <- champ.DMP(beta = myCombatPD,
                         pheno = pD$sex,
                         compare.group = NULL,
                         adjPVal = 1,
                         adjust.method = "BH",
                         arraytype = "450K")
head(myDMPPD1_128[[1]]) 
dim(myDMPPD1_128$M_to_F)
save(myDMPPD1_128,file = 'GSE111223_champ_myDMPPD1_128.Rdata')
load(file = 'GSE111223_champ_myDMPPD1_128.Rdata')


library(qqman)
qq(myDMPPD1_128$M_to_F$P.Value)
p_value=myDMPPD1_128$M_to_F$P.Value
z = qnorm(p_value/ 2)
lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)
lambda = round(median(z^2, na.rm = TRUE) / qchisq(0.5, 1), 3)
lambda = round(median(z^2, na.rm = TRUE) / 0.456, 3)



## data$tvalue <- data$Estimate / data$StdErr
## data$zvalue <- qnorm(pt(data$tvalue, df))
## data$chisq <- (data$zvalue) ^ 2
Estimate = myDMPPD1_128$M_to_F$logFC
StdErr = Estimate/myDMPPD1_128$M_to_F$t
zvalue = qnorm(myDMPPD1_128$M_to_F$P.Value / 2)
inflationFactor = median(zvalue^2,na.rm = TRUE) / qchisq(0.5, 1)
print("lambda")
print(inflationFactor)

# genome-wide sig cpgs
sig <- ifelse(myDMPPD1_128$M_to_F$P.Value < 2.4e-7, 1, 0)
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
myDMPPD1_128_bacon <- data.frame(
  myDMPPD1_128,
  Estimate.bacon = bacon::es(bc),
  StdErr.bacon = bacon::se(bc),
  pValue.bacon = pval(bc),
  fdr.bacon = p.adjust(pval(bc), method = "fdr"),
  stringsAsFactors = FALSE)

class(myDMPPD1_128_bacon)  # "data.frame"
keep <- myDMPPD1_128_bacon$M_to_F.adj.P.Val < 0.05  
myDMPPD128_0.05_bacon = myDMPPD1_128_bacon[keep,]
write.table(myDMPPD128_0.05_bacon, file = "myDMPPD128_0.05_bacon.csv", sep = ",")





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

myLoadHC <-champ.load("GSE111223/idat",
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



champ.QC(beta = myLoadHC$beta,
         pheno=myLoadHC$pd$sex,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir="./CHAMP_QCimages_131HC/")

myNormHC <- champ.norm(beta=myLoadHC$beta,
                       rgSet=myLoadHC$rgSet,
                       mset=myLoadHC$mset,
                       resultsDir="./CHAMP_Normalization_131HC/",
                       method="BMIQ",
                       plotBMIQ=FALSE,
                       arraytype="450K",
                       cores=3)
class(myNormHC)   #"matrix" "array" 
dim(myNormHC)  #[1]384201    131
pD=myLoadHC$pd
pD$Slide<-as.character(as.numeric(pD$Slide))
pD$age<-as.numeric(as.character(pD$age))
save(myNormHC,pD,file = 'GSE111223_champ_myNormHC131.Rdata')  
load(file = 'GSE111223_champ_myNormHC131.Rdata')

champ.SVD(beta = myNormHC %>% as.data.frame(),
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_131HC/")


?champ.runCombat
myCombatHC <- champ.runCombat(beta=myNormHC,
                              pd=pD,
                              variablename="sex",
                              batchname=c("Slide"),
                              logitTrans=TRUE)

champ.SVD(beta = myCombatHC %>% as.data.frame(),
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_131HC_AftercombatSlide/")

myCombatHC <- champ.runCombat(beta=myCombatHC,
                              pd=pD,
                              variablename="sex",
                              batchname=c("Array"),
                              logitTrans=TRUE)
save(myCombatPD,pD,file = 'GSE145361_champ_myCombatPD813.Rdata')
load(file = 'GSE145361_champ_myCombatPD813.Rdata')


champ.SVD(beta = myCombatHC %>% as.data.frame(),
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_131HC_AftercombatSlideArray/")


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
myDMPHC0.05_131 <- champ.DMP(beta = myCombatHC,
                     pheno = pD$sex,
                     compare.group = NULL,
                     adjPVal = 0.05,
                     adjust.method = "BH",
                     arraytype = "450K")
head(myDMPHC0.05_131[[1]]) 
dim(myDMPHC0.05_131$M_to_F)
save(myDMPHC0.05_131,file = 'GSE111223_champ_myDMPHC0.05_131.Rdata')
load(file = 'GSE111223_champ_myDMPHC0.05_131.Rdata')
write.table(myDMPHC0.05_131, file = "myDMPHC0.05_131.csv", sep = ",")

myDMP <- champ.DMP(beta = myNormHC,
                   pheno = pD$sex,
                   compare.group = NULL,
                   adjPVal = 0.05,
                   adjust.method = "BH",
                   arraytype = "450K")
head(myDMP[[1]])
dim(myDMP$M_to_F)
rm(myDMP)


myDMPHC1_131 <- champ.DMP(beta = myCombatHC,
                         pheno = pD$sex,
                         compare.group = NULL,
                         adjPVal = 1,
                         adjust.method = "BH",
                         arraytype = "450K")
head(myDMPHC1_131[[1]]) 
dim(myDMPHC1_131$M_to_F)
save(myDMPHC1_131,file = 'GSE111223_champ_myDMPHC1_131.Rdata')
load(file = 'GSE111223_champ_myDMPHC1_131.Rdata')


library(qqman)
qq(myDMPHC1_131$M_to_F$P.Value)
p_value=myDMPHC1_131$M_to_F$P.Value
z = qnorm(p_value/ 2)
lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)
lambda = round(median(z^2, na.rm = TRUE) / qchisq(0.5, 1), 3)
lambda = round(median(z^2, na.rm = TRUE) / 0.456, 3)



## data$tvalue <- data$Estimate / data$StdErr
## data$zvalue <- qnorm(pt(data$tvalue, df))
## data$chisq <- (data$zvalue) ^ 2
Estimate = myDMPHC1_131$M_to_F$logFC
StdErr = Estimate/myDMPHC1_131$M_to_F$t
zvalue = qnorm(myDMPHC1_131$M_to_F$P.Value / 2)
inflationFactor = median(zvalue^2,na.rm = TRUE) / qchisq(0.5, 1)
print("lambda")
print(inflationFactor)

# genome-wide sig cpgs
sig <- ifelse(myDMPHC1_131$M_to_F$P.Value < 2.4e-7, 1, 0)
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
myDMPHC1_131_bacon <- data.frame(
  myDMPHC1_131,
  Estimate.bacon = bacon::es(bc),
  StdErr.bacon = bacon::se(bc),
  pValue.bacon = pval(bc),
  fdr.bacon = p.adjust(pval(bc), method = "fdr"),
  stringsAsFactors = FALSE)

class(myDMPHC1_131_bacon)  # "data.frame"
keep <- myDMPHC1_131_bacon$M_to_F.adj.P.Val < 0.05  
myDMPHC131_0.05_bacon = myDMPHC1_131_bacon[keep,]
write.table(myDMPHC131_0.05_bacon, file = "myDMPHC131_0.05_bacon.csv", sep = ",")




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




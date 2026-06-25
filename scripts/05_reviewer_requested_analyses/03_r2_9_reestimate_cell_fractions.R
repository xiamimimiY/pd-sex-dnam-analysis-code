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



install.packages("vctrs", version = "0.5.2", repos = "http://cran.us.r-project.org")

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

targetsall <- read.metharray.sheet("ppmi_140_idat all", pattern="ppmi_140_samplesheet_sPD_BL.csv") 
targetsall <- targetsall[-93, ]


rgSetall <- read.metharray.exp("ppmi_140_idat all",targets=targetsall, force=TRUE)
dim(rgSetall) #[1]  1051539     214
class(rgSetall)  #[1] "RGChannelSet" attr(,"package")  [1] "minfi"
save(rgSetall,file = 'ppmi_140_minfi_rgSetall_sPD_BL.Rdata') 
load(file = 'ppmi_140_minfi_rgSetall_sPD_BL.Rdata')


detP <- detectionP(rgSetall)
dim(detP)  #[1] 865859    214
head(detP)
?ColorBrewer
library(RColorBrewer)
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targetsall$PATNO)], las=2, 
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
dim(Msetall)  #[1] 865859    214
GMsetall <- mapToGenome(Msetall)
dim(GMsetall)  #[1]865859    214

PreSex <- getSex(object = GMsetall, cutoff = -2)
class(PreSex) #DFrame

targetsall$Sex == PreSex$predictedSex
all(targetsall$Sex == PreSex$predictedSex)
PreSex$predictedSex == ppmi_140_samplesheet_sPD_V08$Sex
all(PreSex$predictedSex == ppmi_140_samplesheet_sPD_V08$Sex)
rm (ppmi_140_samplesheet_sPD_V08)
targetsall <- targetsall[-47, ]




library(FlowSorted.Blood.450k)
cellCounts <- estimateCellCounts(rgSetall)
save(cellCounts,file = 'GSE111629_minfi_cellCounts.csv')



myLoadPD <-champ.load("ppmi_140_idat all",
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
           arraytype="EPIC")
dim(myLoadPD$beta) #[1] 736133    213
save(myLoadPD,file = 'ppmi140_BL_champ_myLoadPD213.Rdata')
load(file = 'ppmi140_BL_champ_myLoadPD213.Rdata')



champ.QC(beta = myLoadPD$beta,
         pheno=myLoadPD$pd$Sex,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir="./CHAMP_QCimages_BL_PD213/")

myNormPD <- champ.norm(beta=myLoadPD$beta,
                       rgSet=myLoadPD$rgSet,
                       mset=myLoadPD$mset,
                       resultsDir="./CHAMP_Normalization_BL_PD213/",
                       method="BMIQ",
                       plotBMIQ=FALSE,
                       arraytype="EPIC",
                       cores=3)
class(myNormPD)   #"matrix" "array" 
dim(myNormPD)  #[1]736133    213
pD=myLoadPD$pd
pD$Slide<-as.character(as.numeric(pD$Slide))
class(pD$Age) # "numeric"
#save(myNormPD,pD,file = 'GSE145361_champ_myNormPD813.Rdata')   #12：59-


champ.SVD(beta = myNormPD,
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_V4_PD190/")

class(myNormPD) %>% length()   #equals 2
class(myNormPD) == "data.frame"   #So it will cause problems:  the condition has length > 1 and only the first element will be used
## u can code this champ.SVD(beta=myNorm %>% as.data.frame(), pd=myLoad$pd) to avoid it

champ.SVD(beta = myNormPD %>% as.data.frame(),
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_BL_PD213/")


?champ.runCombat
myCombatPD <- champ.runCombat(beta=myNormPD,
                              pd=pD,
                              variablename="Sex",
                              batchname=c("Slide"),
                              logitTrans=TRUE)
dim(myCombatPD) #736133    213

champ.SVD(beta = myCombatPD %>% as.data.frame(),
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_BL_PD213_AftercombatSlide/")

myCombatPD <- champ.runCombat(beta=myCombatPD,
                              pd=pD,
                              variablename="Sex",
                              batchname=c("Array"),
                              logitTrans=TRUE)


champ.SVD(beta = myCombatPD %>% as.data.frame(),
          rgSet=NULL,
          pd=pD,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_BL_PD213_AftercombatSlideArray/")


?champ.refbase
myRefbasePD <- champ.refbase(beta = myCombatPD,
              arraytype="EPIC")
class(myRefbasePD)  #[1] "list"
class(myRefbasePD$CorrectedBeta)   #[1] "matrix" "array" 
class(myRefbasePD$CellFraction)  #[1] "matrix" "array" 
class(myNormPD)    #[1] "matrix" "array"
save(myRefbasePD,pD,file = 'ppmi140_champ_myRefbase_BL_PD213.Rdata')   
load(file = 'ppmi140_champ_myRefbase_BL_PD213.Rdata')


?champ.DMP
myDMPPD0.05_BL_PD213 <- champ.DMP(beta = myRefbasePD$CorrectedBeta,
                     pheno = pD$Sex,
                     compare.group = NULL,
                     adjPVal = 0.05,
                     adjust.method = "BH",
                     arraytype = "EPIC")
head(myDMPPD0.05_BL_PD213[[1]]) 
dim(myDMPPD0.05_BL_PD213$M_to_F)
save(myDMPPD0.05_BL_PD213,file = 'ppmi140_champ_myDMP0.05_BL_PD213.Rdata')
load(file = 'ppmi140_champ_myDMP0.05_V8_PD197.Rdata')
write.table(myDMPPD0.05_BL_PD213, file = "myDMP0.05_BL_PD213.csv", sep = ",")

myDMP <- champ.DMP(beta = myNormPD,
                   pheno = pD$Sex,
                   compare.group = NULL,
                   adjPVal = 0.05,
                   adjust.method = "BH",
                   arraytype = "EPIC")
head(myDMP[[1]])
dim(myDMP$M_to_F)
myDMP <- champ.DMP(beta = myCombatPD,
                   pheno = pD$Sex,
                   compare.group = NULL,
                   adjPVal = 0.05,
                   adjust.method = "BH",
                   arraytype = "EPIC")
head(myDMP[[1]])
dim(myDMP$M_to_F)
rm(myDMP)


myDMP1_BL_PD213 <- champ.DMP(beta = myRefbasePD$CorrectedBeta,
                         pheno = pD$Sex,
                         compare.group = NULL,
                         adjPVal = 1,
                         adjust.method = "BH",
                         arraytype = "EPIC")
head(myDMP1_BL_PD213[[1]]) 
dim(myDMP1_BL_PD213$M_to_F)
save(myDMP1_BL_PD213,file = 'ppmi140_champ_myDMP1_BL_PD213.Rdata')
load(file = 'ppmi140_champ_myDMP1_BL_PD213.Rdata')


library(qqman)
qq(myDMP1_BL_PD213$M_to_F$P.Value)  
p_value=myDMP1_BL_PD213$M_to_F$P.Value
z = qnorm(p_value/ 2)
lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)
lambda = round(median(z^2, na.rm = TRUE) / qchisq(0.5, 1), 3)
lambda = round(median(z^2, na.rm = TRUE) / 0.456, 3)

Estimate = myDMP1_BL_PD213$M_to_F$logFC
StdErr = Estimate/myDMP1_BL_PD213$M_to_F$t
zvalue = qnorm(myDMP1_BL_PD213$M_to_F$P.Value / 2)
inflationFactor = median(zvalue^2,na.rm = TRUE) / qchisq(0.5, 1)
print("lambda")
print(inflationFactor)
rm(zvalue)
# genome-wide sig cpgs
sig <- ifelse(myDMP1_V8_PD197$M_to_F$P.Value < 2.4e-7, 1, 0)
table(sig)
rm(sig)

### 2. bacon analysis
library(bacon)
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
myDMP1_BL_PD213_bacon <- data.frame(
  myDMP1_BL_PD213,
  Estimate.bacon = bacon::es(bc),
  StdErr.bacon = bacon::se(bc),
  pValue.bacon = pval(bc),
  fdr.bacon = p.adjust(pval(bc), method = "fdr"),
  stringsAsFactors = FALSE)

class(myDMP1_BL_PD213_bacon)  # "data.frame"
keep <- myDMP1_BL_PD213_bacon$M_to_F.adj.P.Val < 0.05  
myDMP0.05_BL_PD213_bacon = myDMP1_BL_PD213_bacon[keep,]
dim(myDMP0.05_BL_PD213_bacon) #[1]183643     24
write.table(myDMP0.05_BL_PD213_bacon, file = "myDMP0.05_BL_PD213_bacon.csv", sep = ",")





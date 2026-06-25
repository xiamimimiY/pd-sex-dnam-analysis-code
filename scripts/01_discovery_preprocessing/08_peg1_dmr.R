#2023-5-1 
rm(list = ls())
options(stringsAsFactors = F)


library(ChAMP)
load(file = 'GSE111629_champ_myRefbasePD330.Rdata')
myRefbasePD$CorrectedBeta[1:5, 1:5]

pD <- read.csv("330PDsamplesheet.csv", header = TRUE, row.names = 1)
head(pD)
class(pD)

?champ.DMR
?DMRcate
class(pD$Sex) #[1] "character"
class(myRefbasePD$CorrectedBeta) # "matrix" "array" 
PDbeta = as.data.frame(myRefbasePD$CorrectedBeta)
class(PDbeta)
Sexgroup = as.matrix(pD$Sex)
rm(Sexgroup)

install.packages("installr")
require(installr) #load/install + load installr
updateR()

myDMR_Bumphunter_PD330 = champ.DMR(beta=myRefbasePD$CorrectedBeta,
          pheno=pD$Sex,
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
write.table(myDMR_Bumphunter_PD330, file = "myDMR_Bumphunter_PD330.txt", sep = "\t", row.names = FALSE, quote = FALSE)
read.table("myDMR_Bumphunter_PD330.txt", header = TRUE, sep = "\t")


load(file = 'GSE111629_champ_myNormPD330.Rdata')
class(myNormPD)
class(myRefbasePD$CorrectedBeta)
str(pD$Sex)
rm(PDbeta)


class(myRefbasePD$CorrectedBeta)
min(myRefbasePD$CorrectedBeta)
max(myRefbasePD$CorrectedBeta)

DMR_PD <- champ.DMR(beta=myNormPD,
                    pheno=pD$Sex,
                    compare.group=c("F","M"),
                    arraytype="450K",
                    method = "DMRcate",
                    minProbes=7,
                    adjPvalDmr=0.05,
                    cores=3,
                    rmSNPCH=T,
                    fdr=0.05,
                    dist=2,
                    mafcut=0.05,
                    lambda=1000,
                    C=2)    #18：52-
 
library(DMRcate)
rm(pheno)
?cpg.annotate
library(ExperimentHub)
library(limma)
eh <- ExperimentHub()
FlowSorted.Blood.EPIC <- eh[["EH1136"]]
tcell <- FlowSorted.Blood.EPIC[,colData(FlowSorted.Blood.EPIC)$CD4T==100 |
                                 colData(FlowSorted.Blood.EPIC)$CD8T==100]
head(tcell)
detP <- minfi::detectionP(tcell)
remove <- apply(detP, 1, function (x) any(x > 0.01))
tcell <- tcell[!remove,]
tcell <- minfi::preprocessFunnorm(tcell)
#Subset to chr2 only
tcell <- tcell[seqnames(tcell) == "chr2",]
tcellms <- minfi::getM(tcell)
tcellms.noSNPs <- rmSNPandCH(tcellms, dist=2, mafcut=0.05)
tcell$Replicate[tcell$Replicate==""] <- tcell$Sample_Name[tcell$Replicate==""]
tcellms.noSNPs <- avearrays(tcellms.noSNPs, tcell$Replicate)
tcell <- tcell[,!duplicated(tcell$Replicate)]
tcell <- tcell[rownames(tcellms.noSNPs),]
colnames(tcellms.noSNPs) <- colnames(tcell)
assays(tcell)[["M"]] <- tcellms.noSNPs
assays(tcell)[["Beta"]] <- minfi::ilogit2(tcellms.noSNPs)
type <- factor(tcell$CellType)
design <- model.matrix(~type) 
myannotation <- cpg.annotate("array", tcell, arraytype = "EPIC",
                             analysis.type="differential", design=design, coef=2)
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
groups <- c(CD8T="magenta", CD4T="forestgreen")
cols <- groups[as.character(type)]
DMR.plot(ranges=results.ranges, dmr=1, CpGs=minfi::getBeta(tcell), what="Beta",
         arraytype = "EPIC", phen.col=cols, genome="hg19")


DMRtcell = dmrcate(myannotation,
        lambda = 1000,
        C=NULL,
        pcutoff = "fdr",
        consec = FALSE,
        conseclambda = 10,
        betacutoff = NULL,
        min.cpgs = 2)
head(tcell)
#————————————————————————————————————————————————————————————————————————————
#————————————————————————————————————————————————————————————————————————————
load(file = 'GSE111629_champ_myRefbasePD330.Rdata')
myRefbasePD$CorrectedBeta[1:5, 1:5]

pD <- read.csv("330PDsamplesheet.csv", header = TRUE, row.names = 1)
head(pD)
class(pD)
rm()

sextype <- factor(pD$Sex)
phenoPDSex  <- model.matrix(~sextype) 

myannotation_PD <- cpg.annotate("array", myRefbasePD$CorrectedBeta, what="Beta", arraytype = "450K",
                                analysis.type="differential", design=phenoPDSex,fdr = 0.05, coef=2)

DMR_PD = dmrcate(myannotation_PD,
                 lambda = 1000,
                 C=NULL,
                 pcutoff = "fdr",
                 consec = FALSE,
                 conseclambda = 10,
                 betacutoff = NULL,
                 min.cpgs = 7) # DMResults object with 742 DMRs.
class(DMR_PD)
?extractRanges
myDMR_PD = extractRanges(DMR_PD, genome = "hg19") 
class(myDMR_PD) #[1] "GRanges" attr(,"package") [1] "GenomicRanges"
save(myDMR_PD,pD,file = 'GSE111629_DMRcate_PD330.Rdata')   
load(file = 'GSE111629_DMRcate_PD330.Rdata')
write.table(myDMR_PD, file = "myDMR_DMRcate_PD330.txt", sep = "\t", row.names = FALSE, quote = FALSE)
myDMR_PD <- read.table("myDMR_DMRcate_PD330.txt", header = TRUE, sep = "\t")

groups <- c(M="Male", F="Female")
cols <- groups[as.character(type)]
DMR.plot(ranges=results.ranges, dmr=1, CpGs=minfi::getBeta(tcell), what="Beta",
         arraytype = "EPIC", phen.col=cols, genome="hg19")

library(missMethyl)
enrichment_GO <- goregion(myDMR_PD[1:100], all.cpg = rownames(myRefbasePD$CorrectedBeta),
                          collection = "GO", array.type = "450K")
enrichment_GO <- enrichment_GO[order(enrichment_GO$P.DE),]
head(as.matrix(enrichment_GO), 10)
write.table(enrichment_GO, file = "GO_DMRcate_PD330.txt", sep = "\t", row.names = TRUE, quote = FALSE)


#——————————————————————————————————————————————————————————————————————————————————————
rm(list = ls())
options(stringsAsFactors = F)


library(ChAMP)
load(file = 'GSE111629_champ_myRefbaseHC236.Rdata')
myRefbaseHC$CorrectedBeta[1:5, 1:5]

myDMR_Bumphunter_HC236 = champ.DMR(beta=myRefbaseHC$CorrectedBeta,
                                   pheno=pD$Sex,
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
write.table(myDMR_Bumphunter_HC236, file = "myDMR_Bumphunter_HC236.txt", sep = "\t", row.names = FALSE, quote = FALSE)
read.table("myDMR_Bumphunter_HC236.txt", header = TRUE, sep = "\t")



sextype <- factor(pD$Sex)
phenoPDSex  <- model.matrix(~sextype) 

myannotation_HC <- cpg.annotate("array", myRefbaseHC$CorrectedBeta, what="Beta", arraytype = "450K",
                                analysis.type="differential", design=phenoPDSex,fdr = 0.05, coef=2)

DMR_HC = dmrcate(myannotation_HC,
                 lambda = 1000,
                 C=NULL,
                 pcutoff = "fdr",
                 consec = FALSE,
                 conseclambda = 10,
                 betacutoff = NULL,
                 min.cpgs = 7) # DMResults object with 332 DMRs.
class(DMR_HC)
?extractRanges
myDMR_HC = extractRanges(DMR_HC, genome = "hg19") 
class(myDMR_HC) #[1] "GRanges" attr(,"package") [1] "GenomicRanges"
save(myDMR_HC,pD,file = 'GSE111629_DMRcate_HC236.Rdata')   
load(file = 'GSE111629_DMRcate_HC236.Rdata')
write.table(myDMR_HC, file = "myDMR_DMRcate_HC236.txt", sep = "\t", row.names = FALSE, quote = FALSE)
myDMR_HC <- read.table("myDMR_DMRcate_HC236.txt", header = TRUE, sep = "\t")

groups <- c(M="Male", F="Female")
cols <- groups[as.character(type)]
DMR.plot(ranges=results.ranges, dmr=1, CpGs=minfi::getBeta(tcell), what="Beta",
         arraytype = "EPIC", phen.col=cols, genome="hg19")

library(missMethyl)
enrichment_GO <- goregion(myDMR_HC[1:100], all.cpg = rownames(myRefbaseHC$CorrectedBeta),
                          collection = "GO", array.type = "450K")
enrichment_GO <- enrichment_GO[order(enrichment_GO$P.DE),]
head(as.matrix(enrichment_GO), 10)
write.table(enrichment_GO, file = "GO_DMRcate_HC236.txt", sep = "\t", row.names = TRUE, quote = FALSE)






###2023-5-26
library(DMRcate)
load(file = 'GSE111629_champ_myRefbasePD330.Rdata')

pD <- read.csv("330PDsamplesheet.csv", header = TRUE, row.names = 1)

sextype <- factor(pD$Sex)
phenoPDSex  <- model.matrix(~sextype) 

myannotation_PD <- cpg.annotate("array", myRefbasePD$CorrectedBeta, what="Beta", arraytype = "450K",
                                analysis.type="differential", design=phenoPDSex,fdr = 0.05, coef=2)

DMR_PD = dmrcate(myannotation_PD,
                 lambda = 1000,
                 C=NULL,
                 pcutoff = "fdr",
                 consec = FALSE,
                 conseclambda = 10,
                 betacutoff = NULL,
                 min.cpgs = 4) # DMResults object with 1344 DMRs.

DMR_PD1 = dmrcate(myannotation_PD,
                 lambda = 500,
                 C=5,
                 pcutoff = "fdr",
                 consec = FALSE,
                 conseclambda = 10,
                 betacutoff = NULL,
                 min.cpgs = 4) # DMResults object with 1123 DMRs.
rm(DMR_PD)
DMR_PD = DMR_PD1
rm(DMR_PD1)

myDMR_PD = extractRanges(DMR_PD, genome = "hg19") 
class(myDMR_PD) #[1] "GRanges" attr(,"package") [1] "GenomicRanges"
write.table(myDMR_PD, file = "myDMR_DMRcate_PD330_mincpg4.txt", sep = "\t", row.names = FALSE, quote = FALSE)

myDMR7_PD <- read.table("myDMR_DMRcate_PD330.txt", header = TRUE, sep = "\t")  #742
myDMR4_PD <- read.table("myDMR_DMRcate_PD330_mincpg4.txt", header = TRUE, sep = "\t")  #1123
new_dmr <- myDMR4_PD[myDMR4_PD$no.cpgs >= 7, ] #329 
rm(new_dmr)

#————————————————————————————————————————————————————————————————————————————————
rm(list = ls())
options(stringsAsFactors = F)

load(file = 'GSE111629_champ_myRefbaseHC236.Rdata')

sextype <- factor(pD$Sex)
phenoPDSex  <- model.matrix(~sextype) 

myannotation_HC <- cpg.annotate("array", myRefbaseHC$CorrectedBeta, what="Beta", arraytype = "450K",
                                analysis.type="differential", design=phenoPDSex,fdr = 0.05, coef=2)

DMR_HC = dmrcate(myannotation_HC,
                 lambda = 1000,
                 C=NULL,
                 pcutoff = "fdr",
                 consec = FALSE,
                 conseclambda = 10,
                 betacutoff = NULL,
                 min.cpgs = 4)

DMR_HC = dmrcate(myannotation_HC,
                  lambda = 500,
                  C=5,
                  pcutoff = "fdr",
                  consec = FALSE,
                  conseclambda = 10,
                  betacutoff = NULL,
                  min.cpgs = 4) # DMResults object with 465 DMRs.
myDMR_HC = extractRanges(DMR_HC, genome = "hg19") 
write.table(myDMR_HC, file = "myDMR_DMRcate_HC236_mincpg4.txt", sep = "\t", row.names = FALSE, quote = FALSE)
myDMR_HC <- read.table("myDMR_DMRcate_HC236_mincpg4.txt", header = TRUE, sep = "\t")

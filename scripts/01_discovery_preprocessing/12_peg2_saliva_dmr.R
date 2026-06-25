
rm(list = ls())
options(stringsAsFactors = F)

library(ChAMP)
class(pD) # "DFrame"
pD <- as.data.frame(pD)

myDMR_Bumphunter_PD128 = champ.DMR(beta=myCombatPD,
                                   pheno=pD$sex,
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

write.table(myDMR_Bumphunter_PD128, file = "myDMR_Bumphunter_salivaPD128.txt", sep = "\t", row.names = FALSE, quote = FALSE)
read.table("myDMR_Bumphunter_salivaPD128.txt", header = TRUE, sep = "\t")


library(DMRcate)
load(file = 'GSE145361_champ_myCombatPD813.Rdata')
myCombatPD[1:5, 1:5]

pD <- read.csv("330PDsamplesheet.csv", header = TRUE, row.names = 1)
head(pD)
class(pD)

sextype <- factor(pD$sex)
phenoPDSex  <- model.matrix(~sextype) 

myannotation_PD <- cpg.annotate("array", myCombatPD, what="Beta", arraytype = "450K",
                                analysis.type="differential", design=phenoPDSex,fdr = 0.05, coef=2)

DMR_PD = dmrcate(myannotation_PD,
                 lambda = 1000,
                 C=NULL,
                 pcutoff = "fdr",
                 consec = FALSE,
                 conseclambda = 10,
                 betacutoff = NULL,
                 min.cpgs = 7) # DMResults object with 168 DMRs.
class(DMR_PD)  #[1] "DMResults" attr(,"package") [1] "DMRcate"
myDMR_PD = extractRanges(DMR_PD, genome = "hg19") 
class(myDMR_PD) #[1] "GRanges" attr(,"package") [1] "GenomicRanges"
save(myDMR_PD,pD,file = 'GSE111223_DMRcate_PD128.Rdata')   
load(file = 'GSE111223_DMRcate_PD128.Rdata')
write.table(myDMR_PD, file = "myDMR_DMRcate_salivaPD128.txt", sep = "\t", row.names = FALSE, quote = FALSE)
myDMR_PD <- read.table("myDMR_DMRcate_salivaPD128.txt", header = TRUE, sep = "\t")


groups <- c(M="Male", F="Female")
cols <- groups[as.character(sextype)]
DMR.plot(ranges=myDMR_PD, dmr=1, CpGs=myCombatPD, what="Beta",
         arraytype = "450K", phen.col=cols, genome="hg19")

library(missMethyl)
enrichment_GO <- goregion(myDMR_PD[1:100], all.cpg = rownames(myCombatPD),
                          collection = "GO", array.type = "450K")
enrichment_GO <- enrichment_GO[order(enrichment_GO$P.DE),]
head(as.matrix(enrichment_GO), 10)
write.table(enrichment_GO, file = "GO_DMRcate_PD813.txt", sep = "\t", row.names = TRUE, quote = FALSE)



#——————————————————————————————————————————————————————————————————————————————————————
rm(list = ls())
options(stringsAsFactors = F)


library(ChAMP)
load(file = 'GSE145361_champ_myRefbaseHC777.Rdata')
myRefbaseHC[1:5, 1:5]

myDMR_Bumphunter_HC777 = champ.DMR(beta=myRefbaseHC,
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
write.table(myDMR_Bumphunter_HC777, file = "myDMR_Bumphunter_HC777.txt", sep = "\t", row.names = FALSE, quote = FALSE)
read.table("myDMR_Bumphunter_HC777.txt", header = TRUE, sep = "\t")



sextype <- factor(pD$Sex)
phenoPDSex  <- model.matrix(~sextype) 

myannotation_HC <- cpg.annotate("array", myRefbaseHC, what="Beta", arraytype = "450K",
                                analysis.type="differential", design=phenoPDSex,fdr = 0.05, coef=2)

DMR_HC = dmrcate(myannotation_HC,
                 lambda = 1000,
                 C=NULL,
                 pcutoff = "fdr",
                 consec = FALSE,
                 conseclambda = 10,
                 betacutoff = NULL,
                 min.cpgs = 7) # DMResults object with 2534 DMRs.
class(DMR_HC)
?extractRanges
myDMR_HC = extractRanges(DMR_HC, genome = "hg19") 
class(myDMR_HC) #[1] "GRanges" attr(,"package") [1] "GenomicRanges"
save(myDMR_HC,pD,file = 'GSE145361_DMRcate_HC777.Rdata')   
load(file = 'GSE145361_DMRcate_HC777.Rdata')
write.table(myDMR_HC, file = "myDMR_DMRcate_HC777.txt", sep = "\t", row.names = FALSE, quote = FALSE)
myDMR_HC <- read.table("myDMR_DMRcate_HC777.txt", header = TRUE, sep = "\t")

groups <- c(M="Male", F="Female")
cols <- groups[as.character(type)]
DMR.plot(ranges=results.ranges, dmr=1, CpGs=minfi::getBeta(tcell), what="Beta",
         arraytype = "EPIC", phen.col=cols, genome="hg19")

library(missMethyl)
enrichment_GO <- goregion(myDMR_HC[1:100], all.cpg = rownames(myRefbaseHC),
                          collection = "GO", array.type = "450K")
enrichment_GO <- enrichment_GO[order(enrichment_GO$P.DE),]
head(as.matrix(enrichment_GO), 10)
write.table(enrichment_GO, file = "GO_DMRcate_HC777.txt", sep = "\t", row.names = TRUE, quote = FALSE)




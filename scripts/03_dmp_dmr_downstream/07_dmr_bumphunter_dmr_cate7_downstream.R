
rm(list = ls())
options(stringsAsFactors = F)

myDMR_Bumphunter_PD813 = read.table("myDMR_Bumphunter_PD813.txt", header = TRUE, sep = "\t")
myDMR_Bumphunter_HC777 = read.table("myDMR_Bumphunter_HC777.txt", header = TRUE, sep = "\t")
myDMR_DMRcate_PD813 = read.table("myDMR_DMRcate_PD813.txt", header = TRUE, sep = "\t")
myDMR_DMRcate_HC777 = read.table("myDMR_DMRcate_HC777.txt", header = TRUE, sep = "\t")
myDMR_Bumphunter_PD330 = read.table("myDMR_Bumphunter_PD330.txt", header = TRUE, sep = "\t")
myDMR_Bumphunter_HC236 = read.table("myDMR_Bumphunter_HC236.txt", header = TRUE, sep = "\t")
myDMR_DMRcate_PD330 = read.table("myDMR_DMRcate_PD330.txt", header = TRUE, sep = "\t")
myDMR_DMRcate_HC236 = read.table("myDMR_DMRcate_HC236.txt", header = TRUE, sep = "\t")

GO_DMRcate_PD813 = read.table("GO_DMRcate_PD813.txt", header = TRUE, sep = "\t")
GO_DMRcate_HC777 = read.table("GO_DMRcate_HC777.txt", header = TRUE, sep = "\t")
GO_DMRcate_PD330 = read.table("GO_DMRcate_PD330.txt", header = TRUE, sep = "\t")
GO_DMRcate_HC236 = read.table("GO_DMRcate_HC236.txt", header = TRUE, sep = "\t")

#_____________________________________________________________________________________
#_____________________________________________________________________________________
myDMR_DMRcate_PD330[1:5,1:5]
cmPD_DMRcate <- merge(myDMR_DMRcate_PD330, myDMR_DMRcate_PD813, by = c("seqnames", "start", "end"))
cmHC_DRMcate <- merge(myDMR_DMRcate_HC236, myDMR_DMRcate_HC777, by = c("seqnames", "start", "end"))
cmPDHC_DMRcate <- merge(cmPD_DMRcate, cmHC_DRMcate, by = c("seqnames", "start", "end"))
library(dplyr)
uniquePD_DMRcate <- anti_join(cmPD_DMRcate, cmPDHC_DMRcate, by = c("seqnames", "start", "end"))

rm(enrichment_GO1)
library(missMethyl)
library(GenomicRanges)
?goregion()
uniquePD_DMRcate_forGO <- makeGRangesFromDataFrame(uniquePD_DMRcate, keep.extra.columns=TRUE)
enrichment_GO <- goregion(uniquePD_DMRcate_forGO, all.cpg = NULL,
                          collection = "GO", array.type = "450K",sig.genes = TRUE)
dim(enrichment_GO) #[1] 22799     6

enrichment_GO <- enrichment_GO[order(enrichment_GO$P.DE),]
enrichment_GO <- enrichment_GO[enrichment_GO$P.DE < 0.05, ]
head(as.matrix(enrichment_GO), 10)
write.table(enrichment_GO, file = "GO_DMRcate_uniquePD.txt", sep = "\t", row.names = TRUE, quote = FALSE)

#_____________________________________________________________________________________
#_____________________________________________________________________________________
colnames(myDMR_Bumphunter_PD813) <- gsub("BumphunterDMR\\.", "", colnames(myDMR_Bumphunter_PD813))
colnames(myDMR_Bumphunter_HC777) <- gsub("BumphunterDMR\\.", "", colnames(myDMR_Bumphunter_HC777))
colnames(myDMR_Bumphunter_PD330) <- gsub("BumphunterDMR\\.", "", colnames(myDMR_Bumphunter_PD330))
colnames(myDMR_Bumphunter_HC236) <- gsub("BumphunterDMR\\.", "", colnames(myDMR_Bumphunter_HC236))

cmPD_Bumphunter <- merge(myDMR_Bumphunter_PD330, myDMR_Bumphunter_PD813, by = c("seqnames", "start", "end"))
cmHC_Bumphunter <- merge(myDMR_Bumphunter_HC236, myDMR_Bumphunter_HC777, by = c("seqnames", "start", "end"))
cmPDHC_Bumphunter <- merge(cmPD_Bumphunter, cmHC_Bumphunter, by = c("seqnames", "start", "end"))
library(dplyr)
uniquePD_Bumphunter <- anti_join(cmPD_Bumphunter, cmPDHC_Bumphunter, by = c("seqnames", "start", "end"))

uniquePD_Bumphunter_forGO <- makeGRangesFromDataFrame(uniquePD_Bumphunter, keep.extra.columns=TRUE)
enrichment_GO1 <- goregion(uniquePD_Bumphunter_forGO, all.cpg = NULL,
                          collection = "GO", array.type = "450K",sig.genes = TRUE)
dim(enrichment_GO1) #[1] 22799     7

enrichment_GO1 <- enrichment_GO1[order(enrichment_GO1$P.DE),]
enrichment_GO1 <- enrichment_GO1[enrichment_GO1$P.DE < 0.05, ]
head(as.matrix(enrichment_GO1), 10)
write.table(enrichment_GO1, file = "GO_Bumphunter_uniquePD.txt", sep = "\t", row.names = TRUE, quote = FALSE)

#_____________________________________________________________________________________
#_____________________________________________________________________________________
uniquePD_Bumphunter_DMRcate_DMR <- merge(uniquePD_Bumphunter, uniquePD_DMRcate, by = c("seqnames", "start", "end"))
uniquePD_Bumphunter_DMRcate_DMR_forGO <- makeGRangesFromDataFrame(uniquePD_Bumphunter_DMRcate_DMR, keep.extra.columns=TRUE)
enrichment_GO2 <- goregion(uniquePD_Bumphunter_DMRcate_DMR_forGO, all.cpg = NULL,
                           collection = "GO", array.type = "450K",sig.genes = TRUE)
dim(enrichment_GO2) #[1] 22799     7
class(enrichment_GO2)

enrichment_GO2 <- enrichment_GO2[order(enrichment_GO2$P.DE),]
enrichment_GO2 <- enrichment_GO2[enrichment_GO2$P.DE < 0.05, ]
head(as.matrix(enrichment_GO2), 10)
write.table(enrichment_GO2, file = "GO_Bumphunter_DMRcate_uniquePD.txt", sep = "\t", row.names = TRUE, quote = FALSE)
write.table(uniquePD_Bumphunter_DMRcate_DMR, file = "uniquePD_Bumphunter_DMRcate_DMR.txt", sep = "\t", row.names = TRUE, quote = FALSE)
write.table(uniquePD_DMRcate, file = "uniquePD_DMRcate_DMR.txt", sep = "\t", row.names = TRUE, quote = FALSE)
write.table(uniquePD_Bumphunter, file = "uniquePD_Bumphunter_DMR.txt", sep = "\t", row.names = TRUE, quote = FALSE)

#_____________________________________________________________________________________
#_____________________________________________________________________________________
GO_DMRcate_PD813$GO <- rownames(GO_DMRcate_PD813)
GO_DMRcate_PD330$GO <- rownames(GO_DMRcate_PD330)
GO_DMRcate_HC777$GO <- rownames(GO_DMRcate_HC777)
GO_DMRcate_HC236$GO <- rownames(GO_DMRcate_HC236)
enrichment_GO$GO <- rownames(enrichment_GO)
enrichment_GO1$GO <- rownames(enrichment_GO1)
enrichment_GO2$GO <- rownames(enrichment_GO2)
GO_DMRcate_PD813 <- GO_DMRcate_PD813[GO_DMRcate_PD813$P.DE < 0.05, ]
GO_DMRcate_PD330 <- GO_DMRcate_PD330[GO_DMRcate_PD330$P.DE < 0.05, ]
GO_DMRcate_HC777 <- GO_DMRcate_HC777[GO_DMRcate_HC777$P.DE < 0.05, ]
GO_DMRcate_HC236 <- GO_DMRcate_HC236[GO_DMRcate_HC236$P.DE < 0.05, ]
cmPD_GO_DMRcate <- merge(GO_DMRcate_PD330, GO_DMRcate_PD813, by = c("GO"))
cmHC_GO_DMRcate <- merge(GO_DMRcate_HC236, GO_DMRcate_HC777, by = c("GO"))
cmPDHC_GO_DMRcate <- merge(cmPD_GO_DMRcate, cmHC_GO_DMRcate, by = c("GO"))
uniquePD_GO_DMRcate <- anti_join(cmPD_GO_DMRcate, cmPDHC_GO_DMRcate, by = c("GO"))

GO <- merge(uniquePD_GO_DMRcate, enrichment_GO, by = c("GO"))
GO <- merge(uniquePD_GO_DMRcate, enrichment_GO1, by = c("GO"))
rm(GO)
#————————————————————————————————————————————————————————————————————————————————
#———————————————————————————————————————————————————————————————————————————————



#_____________________________________________________________________________________
#______________________________________________________________________________________
DMP_uniquePD <- read.csv("uniquePDann0.05.csv")
colnames(DMP_uniquePD)[colnames(DMP_uniquePD) == "F_to_M.gene"] <- "gene"
colnames(uniquePD_DMRcate)[colnames(uniquePD_DMRcate) == "overlapping.genes.x"] <- "gene"
cmgene_DMP_DMR_DMRcate <- merge(DMP_uniquePD, uniquePD_DMRcate, by = c("gene"))
library(GenomicRanges)
library(IRanges)
# convert DMR results to GRanges object
DMR.gr_DMRcate <- GRanges(
  seqnames=uniquePD_DMRcate$seqnames, 
  ranges = IRanges(start=uniquePD_DMRcate$start, end=uniquePD_DMRcate$end)
  )
class(DMR.gr_DMRcate) # "GRanges" attr(,"package") [1] "GenomicRanges"
DMP_uniquePD$F_to_M.CHR <- paste0("chr", DMP_uniquePD$F_to_M.CHR)
# convert DMP results to GRanges object
DMP.gr <- GRanges(
  seqnames=DMP_uniquePD$F_to_M.CHR, 
  ranges = IRanges(start=DMP_uniquePD$F_to_M.MAPINFO, end=DMP_uniquePD$F_to_M.MAPINFO))

# find overlaps
?findOverlaps
cmposition_DMP_DMR_DMRcate <- findOverlaps(query=DMP.gr, subject=DMR.gr_DMRcate)
cmposition_DMP_DMR_DMRcate <- as.data.frame(cmposition_DMP_DMR_DMRcate)
cmposition_DMP <- DMP.gr[cmposition_DMP_DMR_DMRcate$queryHits,]
cmposition_DMR_DMRcate <- DMR.gr_DMRcate[cmposition_DMP_DMR_DMRcate$subjectHits,]
class(cmposition_DMP) #"GRanges" attr(,"package") [1] "GenomicRanges"
cmposition_DMP <- as.data.frame(cmposition_DMP) #[1] 243   5
cmposition_DMR_DMRcate <- as.data.frame(cmposition_DMR_DMRcate) #[1] 243   5
colnames(cmposition_DMP)[colnames(cmposition_DMP) == "start"] <- "F_to_M.MAPINFO"
cmposition_DMP_ann <- merge(cmposition_DMP,DMP_uniquePD,by="F_to_M.MAPINFO")
cmposition_DMR_DMRcate_ann <- merge(cmposition_DMR_DMRcate,uniquePD_DMRcate,by=c("seqnames", "start", "end"))

#————————————————————————————————————————————————————————————————————————
# convert DMR results to GRanges object
DMR.gr_Bumphunter <- GRanges(
  seqnames=uniquePD_Bumphunter$seqnames, 
  ranges = IRanges(start=uniquePD_Bumphunter$start, end=uniquePD_Bumphunter$end)
)
class(DMR.gr_Bumphunter) # "GRanges" attr(,"package") [1] "GenomicRanges"
# find overlaps
cmposition_DMP_DMR_Bumphunter <- findOverlaps(query=DMP.gr, subject=DMR.gr_Bumphunter)
cmposition_DMP_DMR_Bumphunter <- as.data.frame(cmposition_DMP_DMR_Bumphunter)
cmposition_DMP1 <- DMP.gr[cmposition_DMP_DMR_Bumphunter$queryHits,]
cmposition_DMR_Bumphunter <- DMR.gr_Bumphunter[cmposition_DMP_DMR_Bumphunter$subjectHits,]
cmposition_DMP1 <- as.data.frame(cmposition_DMP1) #[1] 48   5
cmposition_DMR_Bumphunter <- as.data.frame(cmposition_DMR_Bumphunter) #[1] 48   5
colnames(cmposition_DMP1)[colnames(cmposition_DMP1) == "start"] <- "F_to_M.MAPINFO"
cmposition_DMP1_ann <- merge(cmposition_DMP1,DMP_uniquePD,by="F_to_M.MAPINFO")
cmposition_DMR_Bumphunter_ann <- merge(cmposition_DMR_Bumphunter,uniquePD_Bumphunter,by=c("seqnames", "start", "end"))

#————————————————————————————————————————————————————
# convert DMR results to GRanges object
DMR.gr_Bumphunter_DMRcate <- GRanges(
  seqnames=uniquePD_Bumphunter_DMRcate_DMR$seqnames, 
  ranges = IRanges(start=uniquePD_Bumphunter_DMRcate_DMR$start, end=uniquePD_Bumphunter_DMRcate_DMR$end)
)
class(DMR.gr_Bumphunter_DMRcate) # "GRanges" attr(,"package") [1] "GenomicRanges"
# find overlaps
cmposition_DMP_DMR_Bumphunter_DMRcate <- findOverlaps(query=DMP.gr, subject=DMR.gr_Bumphunter_DMRcate)
cmposition_DMP_DMR_Bumphunter_DMRcate <- as.data.frame(cmposition_DMP_DMR_Bumphunter_DMRcate)
cmposition_DMP2 <- DMP.gr[cmposition_DMP_DMR_Bumphunter_DMRcate$queryHits,]
cmposition_DMR_Bumphunter_DMRcate <- DMR.gr_Bumphunter_DMRcate[cmposition_DMP_DMR_Bumphunter_DMRcate$subjectHits,]
cmposition_DMP2 <- as.data.frame(cmposition_DMP2) #[1] 6   5
cmposition_DMR_Bumphunter_DMRcate <- as.data.frame(cmposition_DMR_Bumphunter_DMRcate) #[1] 6   5
colnames(cmposition_DMP2)[colnames(cmposition_DMP2) == "start"] <- "F_to_M.MAPINFO"
cmposition_DMP2_ann <- merge(cmposition_DMP2,DMP_uniquePD,by="F_to_M.MAPINFO")
cmposition_DMR_Bumphunter_DMRcate_ann <- merge(cmposition_DMR_Bumphunter_DMRcate,uniquePD_Bumphunter_DMRcate_DMR,by=c("seqnames", "start", "end"))

#——————————————————————————————————————————————————————
write.table(cmposition_DMP_ann, file = "cmposition_DMP_DMRcate_ann.txt", sep = "\t", row.names = TRUE, quote = FALSE) #243
write.table(cmposition_DMP1_ann, file = "cmposition_DMP_Bumphunter_ann.txt", sep = "\t", row.names = TRUE, quote = FALSE) #48
write.table(cmposition_DMP2_ann, file = "cmposition_DMP_Bumphunter_DMRcate_ann.txt", sep = "\t", row.names = TRUE, quote = FALSE) #6





#_____________________________________________________________________________________
#______________________________________________________________________________________
Go_uniquePD = read.table("Go_uniquePD.txt", header = TRUE, sep = "\t")
Go_uniquePD_down = read.table("Go_uniquePD_down.txt", header = TRUE, sep = "\t")
Go_uniquePD_up = read.table("Go_uniquePD_up.txt", header = TRUE, sep = "\t")



#—————————————————————————————————————————————————————————————————————————————————————
genes_uniquePD_DMRcate7 <- uniquePD_DMRcate$gene #188
class(genes_uniquePD_DMRcate7)
genes_uniquePD_DMRcate7 <- na.omit(genes_uniquePD_DMRcate7)
write.table(genes_uniquePD_DMRcate7, file = "genes_uniquePD_DMRcate7.csv",sep = ",", row.names=FALSE,col.names=TRUE)




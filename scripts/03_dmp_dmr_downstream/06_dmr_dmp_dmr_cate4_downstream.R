
###2023-5-27
myDMR_DMRcate4_PD813 = read.table("myDMR_DMRcate_PD813_mincpg4.txt", header = TRUE, sep = "\t") #2736
myDMR_DMRcate4_HC777 = read.table("myDMR_DMRcate_HC777_mincpg4.txt", header = TRUE, sep = "\t") #3298
myDMR_DMRcate4_PD330 = read.table("myDMR_DMRcate_PD330_mincpg4.txt", header = TRUE, sep = "\t") #1123
myDMR_DMRcate4_HC236 = read.table("myDMR_DMRcate_HC236_mincpg4.txt", header = TRUE, sep = "\t") #465
myDMR_combp = read.csv("DMR_combp.csv", header = TRUE)


#_____________________________________________________________________________________
#_____________________________________________________________________________________
myDMR_DMRcate4_PD330[1:5,1:5]
cmPD_DMRcate4 <- merge(myDMR_DMRcate4_PD330, myDMR_DMRcate4_PD813, by = c("seqnames", "start", "end"))
cmHC_DRMcate4 <- merge(myDMR_DMRcate4_HC236, myDMR_DMRcate4_HC777, by = c("seqnames", "start", "end"))
cmPDHC_DMRcate4 <- merge(cmPD_DMRcate4, cmHC_DRMcate4, by = c("seqnames", "start", "end"))
library(dplyr)
uniquePD_DMRcate4 <- anti_join(cmPD_DMRcate4, cmPDHC_DMRcate4, by = c("seqnames", "start", "end"))
write.csv(uniquePD_DMRcate4, file = "uniquePD_DMRcate4.csv", row.names = FALSE)

library(missMethyl)
library(GenomicRanges)
?goregion()
uniquePD_DMRcate4_forGO <- makeGRangesFromDataFrame(uniquePD_DMRcate4, keep.extra.columns=TRUE)
enrichment_GO_DMRcate4 <- goregion(uniquePD_DMRcate4_forGO, all.cpg = NULL,
                          collection = "GO", array.type = "450K",sig.genes = TRUE)
dim(enrichment_GO_DMRcate4) #[1] 22799     7

enrichment_GO_DMRcate4 <- enrichment_GO_DMRcate4[order(enrichment_GO_DMRcate4$P.DE),]
enrichment_GO_DMRcate4 <- enrichment_GO_DMRcate4[enrichment_GO_DMRcate4$P.DE < 0.05, ]
head(as.matrix(enrichment_GO_DMRcate4), 10)
write.table(enrichment_GO_DMRcate4, file = "GO_DMRcate4_uniquePD.txt", sep = "\t", row.names = TRUE, quote = FALSE)



myDMR_combp_forGO <- makeGRangesFromDataFrame(myDMR_combp, keep.extra.columns=TRUE)
enrichment_GO_combp <- goregion(myDMR_combp_forGO, all.cpg = NULL,
                           collection = "GO", array.type = "450K",sig.genes = TRUE)
dim(enrichment_GO_combp) #[1] 22799     7

enrichment_GO_combp <- enrichment_GO_combp[order(enrichment_GO_combp$P.DE),]
enrichment_GO_combp <- enrichment_GO_combp[enrichment_GO_combp$P.DE < 0.05, ]
head(as.matrix(enrichment_GO_combp), 10)
write.table(enrichment_GO_combp, file = "GO_combp_regionbased.txt", sep = "\t", row.names = TRUE, quote = FALSE)

colnames(myDMR_combp)
colnames(myDMR_combp)[colnames(myDMR_combp) == "chrom"] <- "seqnames"
uniquePD_combp_DMRcate_DMR <- merge(myDMR_combp, uniquePD_DMRcate4, by = c("seqnames", "start", "end"))
cm_combp_DMRcate_DMR <- merge(myDMR_combp, myDMR_DMRcate4_PD813, by = c("seqnames", "start", "end"))
cm_combp_DMRcate_DMR <- merge(myDMR_combp, myDMR_DMRcate4_PD813, by = c("seqnames", "start")) 

genes_uniquePD_DMRcate4 <- uniquePD_DMRcate4$overlapping.genes.x #307
class(genes_uniquePD_DMRcate4) #[1] "character"
genes_uniquePD_DMRcate4 <- na.omit(genes_uniquePD_DMRcate4)
write.table(genes_uniquePD_DMRcate4, file = "genes_uniquePD_DMRcate4.csv",sep = ",", row.names=FALSE,col.names=TRUE)










####——————————————————————————————————————————————————————————————————————
genes_uniquePD_DMRcate4 <- uniquePD_DMRcate4$overlapping.genes.x #307
genes_cmPD_DMRcate4 <- cmPD_DMRcate4$overlapping.genes.x #389
genes_cmHC_DRMcate4 <- cmHC_DRMcate4$overlapping.genes.x #168
genes_cmPDHC_DMRcate4 <- cmPDHC_DMRcate4$overlapping.genes.x.x #82
genes_uniquePD_DMRcate4 <- na.omit(genes_uniquePD_DMRcate4) #224
genes_cmPD_DMRcate4 <- na.omit(genes_cmPD_DMRcate4) #289
genes_cmHC_DRMcate4 <- na.omit(genes_cmHC_DRMcate4) #128
genes_cmPDHC_DMRcate4 <- na.omit(genes_cmPDHC_DMRcate4) #65
write.table(genes_cmPD_DMRcate4, file = "genes_cmPD_DMRcate4.csv",sep = ",", row.names=FALSE,col.names=TRUE)
write.table(genes_cmHC_DRMcate4, file = "genes_cmHC_DRMcate4.csv",sep = ",", row.names=FALSE,col.names=TRUE)
write.table(genes_cmPDHC_DMRcate4, file = "genes_cmPDHC_DMRcate4.csv",sep = ",", row.names=FALSE,col.names=TRUE)

library(clusterProfiler)
?enrichGO
Go_uniquePD_DMRcate4 <- enrichGO(gene         = genes_uniquePD_DMRcate4,       # input gene list
                        OrgDb    = "org.Hs.eg.db", 
                        ont         = "ALL",       
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05,      
                        keyType = "SYMBOL")  
print(Go_uniquePD_DMRcate4)

Go_cmPD_DMRcate4 <- enrichGO(gene         = genes_cmPD_DMRcate4,       # input gene list
                                 OrgDb    = "org.Hs.eg.db", # genome background
                                 ont         = "ALL",        # ontology to test (Biological Process)
                                 pAdjustMethod = "BH",      # multiple testing adjustment method
                                 qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                                 keyType = "SYMBOL")  
print(Go_cmPD_DMRcate4)

Go_cmHC_DRMcate4 <- enrichGO(gene         = genes_cmHC_DRMcate4,       # input gene list
                             OrgDb    = "org.Hs.eg.db", # genome background
                             ont         = "ALL",        # ontology to test (Biological Process)
                             pAdjustMethod = "BH",      # multiple testing adjustment method
                             qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                             keyType = "SYMBOL")  
print(Go_cmHC_DRMcate4)

Go_cmPDHC_DMRcate4 <- enrichGO(gene         = genes_cmPDHC_DMRcate4,       # input gene list
                             OrgDb    = "org.Hs.eg.db", # genome background
                             ont         = "ALL",        # ontology to test (Biological Process)
                             pAdjustMethod = "BH",      # multiple testing adjustment method
                             qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                             keyType = "SYMBOL")  
print(Go_cmPDHC_DMRcate4)
write.table(Go_uniquePD_DMRcate4, 'Go_uniquePD_DMRcate4.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(Go_cmPD_DMRcate4, 'Go_cmPD_DMRcate4.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(Go_cmHC_DRMcate4, 'Go_cmHC_DRMcate4.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(Go_cmPDHC_DMRcate4, 'Go_cmPDHC_DMRcate4.txt', sep = '\t', row.names = FALSE, quote = FALSE)


###2023-5-24
dmr <- read.table("out.anno.hg19.bed", sep = "\t", header = TRUE) #1674

dmr1 <- read.table("out1.anno.hg19.bed", sep = "\t", header = TRUE)  #1640

colnames(dmr)

DMR <- dmr[dmr$n_probes > 3, ]  #50
DMR1 <- dmr1[dmr1$n_probes > 3, ]

DMR_combp = write.csv(DMR1, "DMR_combp.csv", row.names = FALSE)

library(clusterProfiler)
genes_combp <- DMR1$refGene_name
DMR_combp_nointergenic <- DMR1[DMR1$refGene_feature != "intergenic", ]
genes_combp_nointergenic <- DMR_combp_nointergenic$refGene_name
Go_combp <- enrichGO(gene         = genes_combp,       # input gene list
                    OrgDb    = "org.Hs.eg.db", # genome background
                    ont         = "ALL",        # ontology to test (Biological Process)
                    pAdjustMethod = "BH",      # multiple testing adjustment method
                    qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                    keyType = "SYMBOL")  
print(Go_combp)
write.table(Go_combp, 'Go_combp.txt', sep = '\t', row.names = FALSE, quote = FALSE)

Go_combp_nointergenic <- enrichGO(gene         = genes_combp_nointergenic,       # input gene list
                     OrgDb    = "org.Hs.eg.db", # genome background
                     ont         = "ALL",        # ontology to test (Biological Process)
                     pAdjustMethod = "BH",      # multiple testing adjustment method
                     qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                     keyType = "SYMBOL")  
print(Go_combp_nointergenic)
write.table(Go_combp_nointergenic, 'Go_combp_nointergenic.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#——————————————————————————————————————————————————————————————————
DMR <- dmr[dmr$n_probes >= 3, ]  #119
DMR3 <- dmr1[dmr1$n_probes >= 3, ]
library(clusterProfiler)
?enrichGO
genes_combp3 <- DMR3$refGene_name
DMR_combp3_nointergenic <- DMR3[DMR3$refGene_feature != "intergenic", ]
genes_combp3_nointergenic <- DMR_combp3_nointergenic$refGene_name
Go_combp3 <- enrichGO(gene         = genes_combp3,       # input gene list
                     OrgDb    = "org.Hs.eg.db", # genome background
                     ont         = "ALL",        # ontology to test (Biological Process)
                     pAdjustMethod = "BH",      # multiple testing adjustment method
                     qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                     keyType = "SYMBOL")  
print(Go_combp3)
write.table(Go_combp3, 'Go_combp3.txt', sep = '\t', row.names = FALSE, quote = FALSE)

Go_combp3_nointergenic <- enrichGO(gene         = genes_combp3_nointergenic,       # input gene list
                                  OrgDb    = "org.Hs.eg.db", # genome background
                                  ont         = "ALL",        # ontology to test (Biological Process)
                                  pAdjustMethod = "BH",      # multiple testing adjustment method
                                  qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                                  keyType = "SYMBOL")  
print(Go_combp3_nointergenic)




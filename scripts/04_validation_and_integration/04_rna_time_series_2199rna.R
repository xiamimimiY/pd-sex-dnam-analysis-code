
####2023-5-31
BL_rna <- read.csv("BL_rna.csv") #55746
V4_rna <- read.csv("V4_rna.csv") #55746
V6_rna <- read.csv("V6_rna.csv") #55746
V8_rna <- read.csv("V8_rna.csv") #55746
uniquePD <- read.csv("uniquePDann0.05.csv", row.names = 1)


colnames(uniquePD)
Anngene <- uniquePD[, c("F_to_M.CHR", "F_to_M.MAPINFO", "F_to_M.gene", "F_to_M.feature", "F_to_M.cgi")]
colnames(Anngene) <- gsub("F_to_M\\.", "", colnames(Anngene))
gene = Anngene$gene
class(gene) #"character" 
gene1 <- gene[gene != ""]
rm(gene)
Anngene <- Anngene[Anngene$gene %in% gene1, ]
duplicate_genes <- Anngene$gene[duplicated(Anngene$gene)]
unique_genes <- Anngene[!duplicated(Anngene$gene), ]


#————————————————————————————————————————————————————————————————————————————————————————————————————
is_x_identical <- identical(BL_rna$X, V4_rna$X) && identical(BL_rna$X, V6_rna$X) && identical(BL_rna$X, V8_rna$X) 
print(is_x_identical)
rm(is_x_identical)

library(biomaRt)
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- BL_rna$X

gene_names <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = ensembl_ids,
                    mart = ensembl)

BL_rna$GeneName <- gene_names$external_gene_name

length(ensembl_ids)
length(unique(ensembl_ids))

length(gene_names$ensembl_gene_id)
length(unique(gene_names$ensembl_gene_id))


#——————————————————————————————————————————————————————————————————————————————————————————————————————
#——————————————————————————————————————————————————————————————————————————————————————————————————————
library(biomaRt)
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
symbol_list <- unique_genes$gene
ensembl_ids <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "external_gene_name",
                     values = symbol_list,
                     mart = ensembl) #1098
ensembl_ids[is.na(ensembl_ids$ensembl_gene_id), "ensembl_gene_id"] <- NA  #1098
ensembl_ids <- aggregate(ensembl_ids["ensembl_gene_id"], by = ensembl_ids["external_gene_name"], FUN = function(x) x[1]) #1098
genes_Tran <- merge(unique_genes, ensembl_ids, by.x = "gene", by.y = "external_gene_name", all.x = TRUE) #1102
write.csv(genes_Tran, file = "genes_Tran.csv", row.names = FALSE)
genes_Tran_add <- read.csv("genes_Tran.csv") 
genes_Tran_add <- na.omit(genes_Tran_add)
colnames(genes_Tran_add) #"gene" "CHR""MAPINFO""feature""cgi""ensembl_gene_id"
colnames(BL_rna)[1] <- "ensembl_gene_id" 
colnames(V4_rna)[1] <- "ensembl_gene_id"
colnames(V6_rna)[1] <- "ensembl_gene_id"
colnames(V8_rna)[1] <- "ensembl_gene_id"
all_rna <- merge(merge(merge(BL_rna, V4_rna, by = "ensembl_gene_id"),
                         V6_rna, by = "ensembl_gene_id"),
                   V8_rna, by = "ensembl_gene_id")
colnames(all_rna)
colnames(all_rna) <- gsub("_rna", "", colnames(all_rna))

TS_rna_male_all = all_rna[, c(1, 2, 4, 6, 8)]
TS_rna_female_all = all_rna[, c(1, 3, 5, 7, 9)]
 
TS_rna_male_ann <- merge(TS_rna_male_all, genes_Tran_add, by = "ensembl_gene_id", all.x = FALSE, all.y = FALSE)
TS_rna_female_ann <- merge(TS_rna_female_all, genes_Tran_add, by = "ensembl_gene_id", all.x = FALSE, all.y = FALSE)
identical(TS_rna_male_ann$ensembl_gene_id, TS_rna_female_ann$ensembl_gene_id) #TRUE
identical(TS_rna_male_ann$gene, TS_rna_female_ann$gene) #TRUE
TS_rna_male = TS_rna_male_ann[, c(1, 2, 3, 4, 5, 6)]
TS_rna_female = TS_rna_female_ann[, c(1, 2, 3, 4, 5)]
rownames(TS_rna_male) <- TS_rna_male[, 1]
TS_rna_male <- TS_rna_male[, -1]
rownames(TS_rna_female) <- TS_rna_female[, 1] 
TS_rna_female <- TS_rna_female[, -1]


#_______________________________________________________________________________
#———————————————————————————————————————————————————————————————————————————————
library(ClusterGVis)
library(Mfuzz)
TS_rna_all <- merge(TS_rna_male, TS_rna_female, by = "row.names", all = TRUE)
class(TS_rna_all) #"data.frame"
rownames(TS_rna_all) <- TS_rna_all[, 6]
TS_rna_all <- TS_rna_all[, -c(1, 6)]
TS_rna_all <- as.matrix(TS_rna_all)

Mset_all <- new("ExpressionSet",exprs = TS_rna_all)
Mset_all <- filter.NA(Mset_all, thres = 0.25)
Mset_all <- fill.NA(Mset_all, mode = 'mean')
Mset_all <- filter.std(Mset_all, min.std = 0)
Mset_all <- standardise(Mset_all)
set.seed(123)
getClusters(exp = na.omit(TS_rna_all))
cluster_Mset_all <- clusterData(exp = Mset_all,
                                cluster.method = "mfuzz",
                                cluster.num = 6) #return a list including wide-shape and long-shape clustered results.
str(cluster_Mset_all)

markGenes = c("CACNG2","CACNA1C","SMAD3","BCL2","MAPK1","RYR2","GRIN2D","HSP90AA1","PPARG","TMEM176B","SKI","NEUROG1","MYF5")

library(org.Hs.eg.db)
enrich <- enrichCluster(object = cluster_Mset_all,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 123)
# check
head(enrich,21)
write.csv(enrich, "enrich.csv", row.names = T)
enrich = read.csv("enrich.csv", header = TRUE, row.names = 1)

pdf('R_all1.pdf',height = 10,width = 11,onefile = F)
visCluster(object = cluster_Mset_all,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           genes.gp = c('italic',fontsize = 9,col = "black"),
           annoTerm.data = enrich,
           line.side = "left",
           go.col = rep(ggsci::pal_d3()(1),each = 1),
           go.size = "pval",
           mulGroup = c(4,4),
           mline.col = c(ggsci::pal_lancet()(3)),
           lgd.label = c("group M", "group F"))
dev.off()


#_______________________________________________________________________________
#———————————————————————————————————————————————————————————————————————————————
str(cluster_Mset_all) #1017 obs. of  11 variables,return a list including wide-shape and long-shape clustered results.
class(cluster_Mset_all) # "list"
cl_rna_all <- as.data.frame(cluster_Mset_all$wide.res)
write.csv(cl_rna_all, "cl_rna_all.csv", row.names = FALSE)


###2023-6-13
# Assuming TS_rna_all is your matrix or array
selected_rows <- TS_rna_all[c("MYF5", "BCL2L10", "CLDN9", "PRSS36", "CDH20", "EVI5L", "PIWIL2", "PXDNL", "TMEM200C", "FAM24B"), ]
selected_rows = as.data.frame(selected_rows)




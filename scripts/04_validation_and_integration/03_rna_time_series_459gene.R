
####2023-6-23
DEG_BL459 = read.csv("DEG_BL459.csv")
colnames(DEG_BL459)[colnames(DEG_BL459) == "ENSEMBL"] <- "ensembl_gene_id"

all_rna <- merge(merge(merge(BL_rna, V4_rna, by = "ensembl_gene_id"),
                         V6_rna, by = "ensembl_gene_id"),
                   V8_rna, by = "ensembl_gene_id")
colnames(all_rna)
colnames(all_rna) <- gsub("_rna", "", colnames(all_rna))

TS_459rna_ann <- merge(all_rna, DEG_BL459, by = "ensembl_gene_id", all.x = FALSE, all.y = FALSE)
colnames(TS_459rna_ann)
row.names(TS_459rna_ann) <- TS_459rna_ann$SYMBOL
TS_459rna <- TS_459rna_ann[, 2:9]
colnames(TS_459rna)
TS_459rna <- TS_459rna[, c("M_BL", "M_V4", "M_V6", "M_V8", "F_BL", "F_V4", "F_V6", "F_V8")]

library(ClusterGVis)
library(Mfuzz)
TS_459rna <- as.matrix(TS_459rna)

Mset_459rna <- new("ExpressionSet",exprs = TS_459rna)
Mset_459rna <- filter.NA(Mset_459rna, thres = 0.25)
Mset_459rna <- fill.NA(Mset_459rna, mode = 'mean')
Mset_459rna <- filter.std(Mset_459rna, min.std = 0)
Mset_459rna <- standardise(Mset_459rna)
set.seed(123)
getClusters(exp = na.omit(TS_459rna))
cluster_Mset_459rna <- clusterData(exp = Mset_459rna,
                                cluster.method = "mfuzz",
                                cluster.num = 6) #return a list including wide-shape and long-shape clustered results.
str(cluster_Mset_459rna)

markGenes_459rna = c("ZNF727","DACT1","DDX43","LTF","CEACAM6","CTSG","OLFM4","MMP8","CEACAM8","DEFA4","OLR1","ABCA13","FIGN","LCNL1","DEFA3","LINC03020")
markGenes_459rna = c("ZNF727","DACT1","DDX43","LTF","CTSG","OLFM4","MMP8","DEFA4")
library(org.Hs.eg.db)
enrich__459rna <- enrichCluster(object = cluster_Mset_459rna,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 123)
# check
head(enrich__459rna)
write.csv(enrich, "enrich.csv", row.names = T)
enrich = read.csv("enrich.csv", header = TRUE, row.names = 1)

pdf('R_459rna1.pdf',height = 10,width = 11,onefile = F)
visCluster(object = cluster_Mset_459rna,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes_459rna,
           markGenes.side = "left",
           genes.gp = c('italic',fontsize = 9,col = "black"),
           annoTerm.data = enrich__459rna,
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




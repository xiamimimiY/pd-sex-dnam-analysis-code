
####2023-5-30
load(file = "ppmi140_champ_myDMP1_BL_PD213.Rdata")
load(file = 'ppmi140_champ_myDMP1_V4_PD190.Rdata')
load(file = 'ppmi140_champ_myDMP1_V6_PD195.Rdata')
load(file = 'ppmi140_champ_myDMP1_V8_PD197.Rdata')
class(myDMP1_V4_PD190) #[1] "list"
myDMP1_BL_PD213 <- do.call(data.frame, myDMP1_BL_PD213)
myDMP1_V4_PD190 <- do.call(data.frame, myDMP1_V4_PD190)
myDMP1_V6_PD195 <- do.call(data.frame, myDMP1_V6_PD195)
myDMP1_V8_PD197 <- do.call(data.frame, myDMP1_V8_PD197)


class(uniquePD0.05) #[1] "character"
BL_2199 <- myDMP1_BL_PD213[uniquePD0.05$MarkerName, ]
V4_2199 <- myDMP1_V4_PD190[uniquePD0.05$MarkerName, ]
V6_2199 <- myDMP1_V6_PD195[uniquePD0.05$MarkerName, ]
V8_2199 <- myDMP1_V8_PD197[uniquePD0.05$MarkerName, ]

rm(myDMP1_BL_PD213)
rm(myDMP1_V4_PD190)
rm(myDMP1_V6_PD195)
rm(myDMP1_V8_PD197)

colnames(BL_2199)
colnames(V4_2199)
colnames(V6_2199)
colnames(V8_2199)
dim(BL_2199)

BL_2199 <- subset(BL_2199, !grepl("NA", rownames(BL_2199))) #2013
V4_2199 <- subset(V4_2199, !grepl("NA", rownames(V4_2199))) #2004
V6_2199 <- subset(V6_2199, !grepl("NA", rownames(V6_2199))) #2012
V8_2199 <- subset(V8_2199, !grepl("NA", rownames(V8_2199))) #2008
common_rows <- Reduce(intersect, list(rownames(BL_2199), rownames(V4_2199), rownames(V6_2199), rownames(V8_2199))) #1985
BL_2199 <- subset(BL_2199, rownames(BL_2199) %in% common_rows)
V4_2199 <- subset(V4_2199, rownames(V4_2199) %in% common_rows)
V6_2199 <- subset(V6_2199, rownames(V6_2199) %in% common_rows)
V8_2199 <- subset(V8_2199, rownames(V8_2199) %in% common_rows)
write.csv(BL_2199, "BL_1985DMP.csv", row.names = T)
write.csv(V4_2199, "V4_1985DMP.csv", row.names = T)
write.csv(V6_2199, "V6_1985DMP.csv", row.names = T)
write.csv(V8_2199, "V8_1985DMP.csv", row.names = T)


TS_cpg_male <- data.frame(
  RowName = rownames(BL_2199),
  M_BL = BL_2199$M_to_F.M_AVG,
  M_V4 = V4_2199$M_to_F.M_AVG,
  M_V6 = V6_2199$M_to_F.M_AVG,
  M_V8 = V8_2199$M_to_F.M_AVG)
TS_cpg_female <- data.frame(
  RowName = rownames(BL_2199),
  F_BL = BL_2199$M_to_F.F_AVG,
  F_V4 = V4_2199$M_to_F.F_AVG,
  F_V6 = V6_2199$M_to_F.F_AVG,
  F_V8 = V8_2199$M_to_F.F_AVG)
class(TS_cpg_male)
head(TS_cpg_female)
rownames(TS_cpg_male) <- TS_cpg_male[, 1]
TS_cpg_male <- TS_cpg_male[, -1]
rownames(TS_cpg_female) <- TS_cpg_female[, 1] 
TS_cpg_female <- TS_cpg_female[, -1]
class(TS_cpg_female)
colnames(TS_cpg_female)

library(Mfuzz)
?mfuzz
TS_cpg_male <- as.matrix(TS_cpg_male)
Mset_male <- new("ExpressionSet",exprs = TS_cpg_male)
boxplot(TS_cpg_male)
Mset_male <- filter.NA(Mset_male, thres = 0.25)
Mset_male <- fill.NA(Mset_male, mode = 'mean')
Mset_male <- filter.std(Mset_male, min.std = 0)
Mset_male <- standardise(Mset_male)
set.seed(123)
cl_Mset_male <- mfuzz(Mset_male, c = 5, m = mestimate(Mset_male))
library(RColorBrewer)
Color <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)
mfuzz.plot(Mset_male,cl_Mset_male,mfrow = c(2,3),new.window = FALSE,time.labels = colnames(TS_cpg_male),colo = Color)
set.seed(123)
library(ClusterGVis)
head(TS_cpg_male)
getClusters(exp = na.omit(TS_cpg_male))
cluster_Mset_male <- clusterData(exp = Mset_male,
                                 cluster.method = "mfuzz",
                                 cluster.num = 6)
str(cluster_Mset_male)

markGenes = rownames(Mset_male)[sample(1:nrow(Mset_male),30,replace = F)]
markGenes = c("cg06376520","cg26998717","cg26385222","cg19318393","cg06756211","cg18310639","cg01963297")
pdf('M_male_1985.pdf',height = 10,width = 6,onefile = F)
visCluster(object = cluster_Mset_male,
           plot.type = "both",
           column_names_rot = 45,
           markGenes = markGenes,
           add.box = T,  
           line.side = "left")
dev.off()

#_______________________________________________________________________________
#_______________________________________________________________________________
TS_cpg_female <- as.matrix(TS_cpg_female)
Mset_female <- new("ExpressionSet",exprs = TS_cpg_female)
Mset_female <- filter.NA(Mset_female, thres = 0.25)
Mset_female <- fill.NA(Mset_female, mode = 'mean')
Mset_female <- filter.std(Mset_female, min.std = 0)
Mset_female <- standardise(Mset_female)
set.seed(123)
getClusters(exp = na.omit(TS_cpg_female))
cluster_Mset_female <- clusterData(exp = Mset_female,
                                   cluster.method = "mfuzz",
                                   cluster.num = 6)
str(cluster_Mset_female)

markGenes = c("cg06376520","cg26998717","cg26385222","cg19318393","cg06756211","cg18310639","cg01963297")
pdf('M_female_1985.pdf',height = 10,width = 6,onefile = F)
visCluster(object = cluster_Mset_female,
           plot.type = "both",
           column_names_rot = 45,
           markGenes = markGenes,
           add.box = T,  
           line.side = "left")
dev.off()

#_______________________________________________________________________________
#_______________________________________________________________________________
TS_cpg_all <- merge(TS_cpg_male, TS_cpg_female, by = "row.names", all = TRUE)
class(TS_cpg_all) #"data.frame"
rownames(TS_cpg_all) <- TS_cpg_all[, 1]
TS_cpg_all <- TS_cpg_all[, -1]
TS_cpg_all <- as.matrix(TS_cpg_all)
write.csv(TS_cpg_all, file = "TS_cpg1985_all.csv", row.names = T) 

Mset_all <- new("ExpressionSet",exprs = TS_cpg_all)
Mset_all <- filter.NA(Mset_all, thres = 0.25)
Mset_all <- fill.NA(Mset_all, mode = 'mean')
Mset_all <- filter.std(Mset_all, min.std = 0)
Mset_all <- standardise(Mset_all)
set.seed(123)
getClusters(exp = na.omit(TS_cpg_all))
cluster_Mset_all <- clusterData(exp = Mset_all,
                                   cluster.method = "mfuzz",
                                   cluster.num = 6)
str(cluster_Mset_all)
write.csv(cluster_Mset_all$wide.res, file = "cluster_Mset_all.csv", row.names = FALSE)


markGenes = c("cg06376520","cg26998717","cg26385222","cg19318393","cg06756211","cg18310639","cg01963297")
markGenes = c("cg06098368","cg03124146","cg12479444","cg20067334","cg14285533","cg26911611","cg26385222","cg06376520","cg26998717","cg17961101")
pdf('M_all1.pdf',height = 10,width = 8,onefile = F)
visCluster(object = cluster_Mset_all,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           genes.gp = c('italic',fontsize = 9,col = "black"),
           line.side = "left",
           mulGroup = c(4,4),
           mline.col = c(ggsci::pal_lancet()(3)),
           lgd.label = c("group M", "group F"),
           set.seed(123))
dev.off()






















#
#————————————————————————————————————————————————————————
class(cluster_Mset_male) #"list"
clus_Mset_male <- do.call(data.frame, cluster_Mset_male$wide.res)
colnames(clus_Mset_male)
colnames(clus_Mset_male)[colnames(clus_Mset_male) == "gene"] <- "MarkerName"
clus_Mset_male <- merge(clus_Mset_male, BL_536[, c("MarkerName", "M_to_F.CHR", "M_to_F.MAPINFO", "M_to_F.gene", "M_to_F.feature", "M_to_F.cgi")], by = "MarkerName", all.x = TRUE)
gene_C1_male <- clus_Mset_male$M_to_F.gene[clus_Mset_male$cluster == 1 & clus_Mset_male$M_to_F.gene != ""]
gene_C2_male <- clus_Mset_male$M_to_F.gene[clus_Mset_male$cluster == 2 & clus_Mset_male$M_to_F.gene != ""] #100
gene_C3_male <- clus_Mset_male$M_to_F.gene[clus_Mset_male$cluster == 3 & clus_Mset_male$M_to_F.gene != ""] #27
gene_C4_male <- clus_Mset_male$M_to_F.gene[clus_Mset_male$cluster == 4 & clus_Mset_male$M_to_F.gene != ""] #90
gene_all_male <- c(gene_C1_male, gene_C2_male, gene_C3_male, gene_C4_male)

library(clusterProfiler)
?enrichGO
GO_all_male <- enrichGO(gene = gene_all_male,   
                        OrgDb    = "org.Hs.eg.db",
                        ont         = "ALL",       
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,      
                        keyType = "SYMBOL")   
print(GO_all_male) #0 enriched terms found
GO_C1_male <- enrichGO(gene = gene_C1_male,  
                      OrgDb    = "org.Hs.eg.db",
                      ont         = "ALL",       
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,      
                     keyType = "SYMBOL")   
print(GO_C1_male) #6 enriched terms found
GO_C1_malea = as.data.frame(GO_C1_male)

GO_C2_male <- enrichGO(gene = gene_C2_male,  
                       OrgDb    = "org.Hs.eg.db",
                       ont         = "ALL",       
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,      
                       keyType = "SYMBOL")   
print(GO_C2_male) #0 enriched terms found


GO_C3_male <- enrichGO(gene = gene_C3_male,  
                       OrgDb    = "org.Hs.eg.db",
                       ont         = "ALL",       
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,      
                       keyType = "SYMBOL")   
print(GO_C3_male) #0 enriched terms found

GO_C4_male <- enrichGO(gene = gene_C4_male,  
                       OrgDb    = "org.Hs.eg.db",
                       ont         = "ALL",       
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,      
                       keyType = "SYMBOL")   
print(GO_C4_male) #0 enriched terms found
#——————————————————————————————————————————————————————————————



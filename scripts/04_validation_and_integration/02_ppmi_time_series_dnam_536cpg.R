
####2023-5-30

library(Mfuzz)

load(file = 'ppmi140_champ_myDMP1_V4_PD190.Rdata')
load(file = 'ppmi140_champ_myDMP1_V6_PD195.Rdata')
load(file = "ppmi140_champ_myDMP1_BL_PD213.Rdata")
load(file = 'ppmi140_champ_myDMP1_V8_PD197.Rdata')
class(myDMP1_V4_PD190) #[1] "list"
myDMP1_V4_PD190 <- do.call(data.frame, myDMP1_V4_PD190)
myDMP1_V6_PD195 <- do.call(data.frame, myDMP1_V6_PD195)
myDMP1_V8_PD197 <- do.call(data.frame, myDMP1_V8_PD197)
myDMP1_BL_PD213 <- do.call(data.frame, myDMP1_BL_PD213)

class(meta_PPMI_BL) #[1] "character"
BL_536 <- PPMI_BL_0.05[meta_PPMI_BL, ]
V4_536 <- myDMP1_V4_PD190[meta_PPMI_BL, ]
V6_536 <- myDMP1_V6_PD195[meta_PPMI_BL, ]
V8_536 <- myDMP1_V8_PD197[meta_PPMI_BL, ]
BL_536_PD213 <- myDMP1_BL_PD213[meta_PPMI_BL, ]
rm(myDMP1_V4_PD190)
rm(myDMP1_V6_PD195)
rm(myDMP1_V8_PD197)
rm(myDMP1_BL_PD213)

colnames(BL_536)
colnames(V4_536)
colnames(V6_536)
colnames(V8_536)
colnames(BL_536_PD213)
TS_cpg_male <- data.frame(
  RowName = rownames(BL_536),
  M_AVG_BL = BL_536$M_to_F.M_AVG,
  M_AVG_V4 = V4_536$M_to_F.M_AVG,
  M_AVG_V6 = V6_536$M_to_F.M_AVG,
  M_AVG_V8 = V8_536$M_to_F.M_AVG)
TS_cpg_female <- data.frame(
  RowName = rownames(BL_536),
  F_AVG_BL = BL_536$M_to_F.F_AVG,
  F_AVG_V4 = V4_536$M_to_F.F_AVG,
  F_AVG_V6 = V6_536$M_to_F.F_AVG,
  F_AVG_V8 = V8_536$M_to_F.F_AVG)
head(TS_cpg_male)
head(TS_cpg_female)
rownames(TS_cpg_male) <- TS_cpg_male[, 1]
TS_cpg_male <- TS_cpg_male[, -1]
rownames(TS_cpg_female) <- TS_cpg_female[, 1] 
TS_cpg_female <- TS_cpg_female[, -1]
class(TS_cpg_female)
colnames(TS_cpg_female)
colnames(TS_cpg_male) <- gsub("_AVG", "", colnames(TS_cpg_male))
colnames(TS_cpg_female) <- gsub("_AVG", "", colnames(TS_cpg_female))


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

#——————————————————————————————————————————————————————————————————————————————————————————————————
#——————————————————————————————————————————————————————————————————————————————————————————————————
library(devtools)
if(! "ClusterGVis" %in% installed.packages()){devtools::install_github("junjunlab/ClusterGVis",dependencies = T)}
# If GitHub installation requires credentials, configure them locally and never commit tokens.
# Configure GitHub credentials locally if needed; do not commit tokens.
library(devtools)
devtools::install_github("junjunlab/ClusterGVis")
library(ClusterGVis)
head(TS_cpg_male)
getClusters(exp = na.omit(TS_cpg_male))
cluster_Mset_male <- clusterData(exp = Mset_male,
                                 cluster.method = "mfuzz",
                                 cluster.num = 5)
str(cluster_Mset_male)
visCluster(object = cluster_Mset_male, plot.type = "line")
visCluster(object = cluster_Mset_male, plot.type = "heatmap")
pdf('M_male.pdf',height = 10,width = 6)
visCluster(object = cluster_Mset_male,
           plot.type = "both",
           column_names_rot = 45)
dev.off()
markGenes = rownames(Mset_male)[sample(1:nrow(Mset_male),30,replace = F)]
markGenes = c("cg06376520","cg26998717","cg26385222","cg19318393","cg06756211","cg18310639","cg01963297")
pdf('M_male_addgene.pdf',height = 10,width = 6,onefile = F)
visCluster(object = cluster_Mset_male,
           plot.type = "both",
           column_names_rot = 45,
           markGenes = markGenes,
           add.box = T,  
           line.side = "left")
dev.off()

#_______________________________________________________________________________
#___________________________________________________________
TS_cpg_female <- as.matrix(TS_cpg_female)
Mset_female <- new("ExpressionSet",exprs = TS_cpg_female)
Mset_female <- filter.NA(Mset_female, thres = 0.25)
Mset_female <- fill.NA(Mset_female, mode = 'mean')
Mset_female <- filter.std(Mset_female, min.std = 0)
Mset_female <- standardise(Mset_female)
set.seed(123)
cl_Mset_female <- mfuzz(Mset_female, c = 5, m = mestimate(Mset_female))
head(TS_cpg_female)
getClusters(exp = na.omit(TS_cpg_female))
cluster_Mset_female <- clusterData(exp = Mset_female,
                                   cluster.method = "mfuzz",
                                   cluster.num = 5)
str(cluster_Mset_female)

markGenes = rownames(Mset_male)[sample(1:nrow(Mset_male),30,replace = F)]
markGenes = c("cg06376520","cg26998717","cg26385222","cg19318393","cg06756211","cg18310639","cg01963297")
pdf('M_female_addgene.pdf',height = 10,width = 6,onefile = F)
visCluster(object = cluster_Mset_female,
           plot.type = "both",
           column_names_rot = 45,
           markGenes = markGenes,
           add.box = T,  
           line.side = "left")
dev.off()








#——————————————————————————————————————————————————————————————————————————————————————
TS_cpg_male_140 <- data.frame(
  RowName = rownames(BL_536_PD213),
  M_BL = BL_536_PD213$M_to_F.M_AVG,
  M_V4 = V4_536$M_to_F.M_AVG,
  M_V6 = V6_536$M_to_F.M_AVG,
  M_V8 = V8_536$M_to_F.M_AVG)
TS_cpg_female_140 <- data.frame(
  RowName = rownames(BL_536_PD213),
  F_BL = BL_536_PD213$M_to_F.F_AVG,
  F_V4 = V4_536$M_to_F.F_AVG,
  F_V6 = V6_536$M_to_F.F_AVG,
  F_V8 = V8_536$M_to_F.F_AVG)
rownames(TS_cpg_male_140) <- TS_cpg_male_140[, 1]
TS_cpg_male_140 <- TS_cpg_male_140[, -1]
rownames(TS_cpg_female_140) <- TS_cpg_female_140[, 1] 
TS_cpg_female_140 <- TS_cpg_female_140[, -1]
class(TS_cpg_female_140)
colnames(TS_cpg_female_140) #"F_BL" "F_V4" "F_V6" "F_V8"

library(Mfuzz)
?mfuzz
TS_cpg_male_140 <- as.matrix(TS_cpg_male_140)
Mset_male_140 <- new("ExpressionSet",exprs = TS_cpg_male_140)
boxplot(TS_cpg_male_140)
Mset_male_140 <- filter.NA(Mset_male_140, thres = 0.25)
Mset_male_140 <- fill.NA(Mset_male_140, mode = 'mean')
Mset_male_140 <- filter.std(Mset_male_140, min.std = 0)
Mset_male_140 <- standardise(Mset_male_140)
set.seed(123)
cl_Mset_male_140 <- mfuzz(Mset_male_140, c = 5, m = mestimate(Mset_male_140))
library(RColorBrewer)
Color <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)
mfuzz.plot(Mset_male_140,cl_Mset_male_140,mfrow = c(2,3),new.window = FALSE,time.labels = colnames(TS_cpg_male_140),colo = Color)


library(ClusterGVis)
head(TS_cpg_male_140)
getClusters(exp = na.omit(TS_cpg_male_140))
cluster_Mset_male_140 <- clusterData(exp = Mset_male_140,
                  cluster.method = "mfuzz",
                  cluster.num = 5)
str(cluster_Mset_male_140)
visCluster(object = cluster_Mset_male_140, plot.type = "line")
visCluster(object = cluster_Mset_male_140, plot.type = "heatmap")
pdf('M_male_140.pdf',height = 10,width = 6)
visCluster(object = cluster_Mset_male_140,
           plot.type = "both",
           column_names_rot = 45)
dev.off()
markGenes = rownames(Mset_male_140)[sample(1:nrow(Mset_male_140),30,replace = F)]
markGenes = c("cg06376520","cg26998717","cg26385222","cg19318393","cg06756211","cg18310639","cg01963297")
pdf('M_male_addgene_140.pdf',height = 10,width = 6,onefile = F)
visCluster(object = cluster_Mset_male_140,
           plot.type = "both",
           column_names_rot = 45,
           markGenes = markGenes,
           add.box = T,  
           line.side = "left")
dev.off()

#_______________________________________________________________________________
#___________________________________________________________
TS_cpg_female_140 <- as.matrix(TS_cpg_female_140)
Mset_female_140 <- new("ExpressionSet",exprs = TS_cpg_female_140)
Mset_female_140 <- filter.NA(Mset_female_140, thres = 0.25)
Mset_female_140 <- fill.NA(Mset_female_140, mode = 'mean')
Mset_female_140 <- filter.std(Mset_female_140, min.std = 0)
Mset_female_140 <- standardise(Mset_female_140)
set.seed(123)
cl_Mset_female <- mfuzz(Mset_female, c = 5, m = mestimate(Mset_female))
head(TS_cpg_female_140)
getClusters(exp = na.omit(TS_cpg_female_140))
cluster_Mset_female_140 <- clusterData(exp = Mset_female_140,
                                 cluster.method = "mfuzz",
                                 cluster.num = 5)
str(cluster_Mset_female_140)

markGenes = c("cg06376520","cg26998717","cg26385222","cg19318393","cg06756211","cg18310639","cg01963297")
pdf('M_female_addgene_140.pdf',height = 10,width = 6,onefile = F)
visCluster(object = cluster_Mset_female_140,
           plot.type = "both",
           column_names_rot = 45,
           markGenes = markGenes,
           add.box = T,  
           line.side = "left")
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



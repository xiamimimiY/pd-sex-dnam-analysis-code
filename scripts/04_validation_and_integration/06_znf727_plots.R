####
rm(list = ls())

gene <- read.csv("path/to/PDSex/eQTMs/eGene_BL.csv", header = TRUE, row.names = 1)
DNAm <- read.csv("path/to/PDSex/eQTMs/DNAm_BL.csv", header = TRUE, row.names = 1)
ZNF727_BL_PDall199 <- gene[rownames(gene) == "ZNF727", ]
ZNF727_BL_PDall199 = as.data.frame(t(ZNF727_BL_PDall199))
row_names <- rownames(ZNF727_BL_PDall199)
new_row_names <- sub("^X", "", row_names)
rownames(ZNF727_BL_PDall199) <- new_row_names

probes6_BL_PDall199 <- DNAm[c("cg06098368", "cg03124146", "cg12479444", "cg20067334", "cg14285533", "cg26911611"), ]
probes6_BL_PDall199 = as.data.frame(t(probes6_BL_PDall199))
row_names <- rownames(probes6_BL_PDall199)
new_row_names <- sub("^X", "", row_names)
rownames(probes6_BL_PDall199) <- new_row_names

write.csv(ZNF727_BL_PDall199, file = "ZNF727gene_BL_PDall199.csv", row.names = T)
write.csv(probes6_BL_PDall199, file = "ZNF727probes6_BL_PDall199.csv", row.names = T)

M_PD_ZNF727 <- read.csv("M_PD_ZNF727.csv", header = TRUE, row.names = 1)
F_PD_ZNF727 <- read.csv("F_PD_ZNF727.csv", header = TRUE, row.names = 1)
M_HC_ZNF727 <- read.csv("M_HC_ZNF727.csv", header = TRUE, row.names = 1)
F_HC_ZNF727 <- read.csv("F_HC_ZNF727.csv", header = TRUE, row.names = 1)

install.packages("ggsci")
library(ggsci)
library(ggplot2)
data <- data.frame(
  Group = factor(c(rep("M_PD", nrow(M_PD_ZNF727)),
                   rep("F_PD", nrow(F_PD_ZNF727)),
                   rep("M_HC", nrow(M_HC_ZNF727)),
                   rep("F_HC", nrow(F_HC_ZNF727))),
                 levels = c("M_PD", "F_PD", "M_HC","F_HC")),
  Expression = c(M_PD_ZNF727$X,
                 F_PD_ZNF727$X,
                 M_HC_ZNF727$X,
                 F_HC_ZNF727$X))


colors <- c("M_PD" = "#4DBBD5B2",
            "F_PD" = "#E64B35B2",
            "M_HC" = "#4DBBD5B2",
            "F_HC" = "#E64B35B2")

plot <- ggplot(data, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values = colors) +
  xlab("Groups") +
  ylab("Expression Level of ZNF727") +
  ggtitle("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())  # Remove minor gridlines
print(plot)
# Save as PDF file
ggsave("ZNF727_exp2.pdf", plot, width = 8, height = 6)


plot <- ggplot(data, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot() +
  scale_fill_jama() +
  xlab("Group") +
  ylab("Normalized counts") +
  ggtitle("Expression Level of ZNF727") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())  # Remove minor gridlines
print(plot)
# Save as PDF file
ggsave("ZNF727_exp4.pdf", plot, width = 8, height = 6)


#——————————————————————————————————————————————————————————————————
library(ggplot2)
library(ggpubr)
library(ggsci)
merged_RNAmDNA <- cbind(ZNF727_BL_PDall199, probes6_BL_PDall199)
c_cg06098368 <- cor(merged_RNAmDNA$ZNF727, merged_RNAmDNA$cg06098368)
p_cg06098368 <- cor.test(merged_RNAmDNA$ZNF727, merged_RNAmDNA$cg06098368)$p.value

plot <- ggscatter(merged_RNAmDNA, x = "cg06098368", y = "ZNF727", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  cor.coef.digits = 2, cor.fontsize = 4,
                  p.value = TRUE, p.value.method = "t.test",
                  p.value.fontsize = 4) +
  xlab("cg06098368 Methylation Level") +
  ylab("ZNF727 Expression Level") +
  scale_fill_jama() +
  ggtitle("")

print(plot)
ggsave("Cor_cg06098368.pdf", plot, width = 8, height = 6)

#cg03124146
plot <- ggscatter(merged_RNAmDNA, x = "cg03124146", y = "ZNF727", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  cor.coef.digits = 2, cor.fontsize = 4,
                  p.value = TRUE, p.value.method = "t.test",
                  p.value.fontsize = 4) +
  xlab("cg03124146 Methylation Level") +
  ylab("ZNF727 Expression Level") +
  scale_fill_jama() +
  ggtitle("")

print(plot)
ggsave("Cor_cg03124146.pdf", plot, width = 8, height = 6)

#cg12479444
plot <- ggscatter(merged_RNAmDNA, x = "cg12479444", y = "ZNF727", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  cor.coef.digits = 2, cor.fontsize = 4,
                  p.value = TRUE, p.value.method = "t.test",
                  p.value.fontsize = 4) +
  xlab("cg12479444 Methylation Level") +
  ylab("ZNF727 Expression Level") +
  scale_fill_jama() +
  ggtitle("")

print(plot)
ggsave("Cor_cg12479444.pdf", plot, width = 8, height = 6)



#cg20067334
plot <- ggscatter(merged_RNAmDNA, x = "cg20067334", y = "ZNF727", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  cor.coef.digits = 2, cor.fontsize = 4,
                  p.value = TRUE, p.value.method = "t.test",
                  p.value.fontsize = 4) +
  xlab("cg20067334 Methylation Level") +
  ylab("ZNF727 Expression Level") +
  scale_fill_jama() +
  ggtitle("")

print(plot)
ggsave("Cor_cg20067334.pdf", plot, width = 8, height = 6)


#cg14285533
plot <- ggscatter(merged_RNAmDNA, x = "cg14285533", y = "ZNF727", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  cor.coef.digits = 2, cor.fontsize = 4,
                  p.value = TRUE, p.value.method = "t.test",
                  p.value.fontsize = 4) +
  xlab("cg14285533 Methylation Level") +
  ylab("ZNF727 Expression Level") +
  scale_fill_jama() +
  ggtitle("")

print(plot)
ggsave("Cor_cg14285533.pdf", plot, width = 8, height = 6)


#cg26911611
plot <- ggscatter(merged_RNAmDNA, x = "cg26911611", y = "ZNF727", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  cor.coef.digits = 2, cor.fontsize = 4,
                  p.value = TRUE, p.value.method = "t.test",
                  p.value.fontsize = 4) +
  xlab("cg26911611 Methylation Level") +
  ylab("ZNF727 Expression Level") +
  scale_fill_jama() +
  ggtitle("")

print(plot)
ggsave("Cor_cg26911611.pdf", plot, width = 8, height = 6)

#————————————————————————————————————————————————————
merged_RNAmDNA$Ave_DNAm <- rowMeans(merged_RNAmDNA[, c("cg06098368", "cg03124146", "cg12479444", "cg20067334", "cg14285533", "cg26911611")])
plot <- ggscatter(merged_RNAmDNA, x = "Ave_DNAm", y = "ZNF727", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  cor.coef.digits = 2, cor.fontsize = 4,
                  p.value = TRUE, p.value.method = "t.test",
                  p.value.fontsize = 4) +
  xlab("Ave_DNAm") +
  ylab("Expression Level of ZNF727") +
  scale_color_jama() +
  ggtitle("")

print(plot)
ggsave("Cor_Ave_DNAm.pdf", plot, width = 8, height = 6)

#__________________________________________________________________________
file_path <- "path/to/PPMI/Me_project 140/ppmi_140_sPD_BL/ppmi_140_clinical_sPD_BL.csv"
sample <- read.csv(file_path, row.names = 1)
merged_RNAmDNA <- merge(merged_RNAmDNA, sample, by = "row.names", all.x = TRUE)
rownames(merged_RNAmDNA) <- merged_RNAmDNA[, 1]
merged_RNAmDNA <- merged_RNAmDNA[, -1]
colnames(merged_RNAmDNA)
merged_RNAmDNA$MDS.UPDRS.part.IA <- as.numeric(as.character(merged_RNAmDNA$MDS.UPDRS.part.IA))
merged_RNAmDNA$MDS.UPDRS.part.IB <- as.numeric(as.character(merged_RNAmDNA$MDS.UPDRS.part.IB))
merged_RNAmDNA$MDS.UPDRS.part.I <- as.numeric(as.character(merged_RNAmDNA$MDS.UPDRS.part.I))
merged_RNAmDNA$MDS.UPDRS.part.II <- as.numeric(as.character(merged_RNAmDNA$MDS.UPDRS.part.II))
merged_RNAmDNA$MDS.UPDRS.part.III <- as.numeric(as.character(merged_RNAmDNA$MDS.UPDRS.part.III))
merged_RNAmDNA$MOCA <- as.numeric(as.character(merged_RNAmDNA$MOCA))
merged_RNAmDNA$GDS <- as.numeric(as.character(merged_RNAmDNA$GDS))
write.csv(merged_RNAmDNA, "merged_RNAmDNA.csv", row.names = TRUE)
library(ggplot2)
library(ggpubr)

# Extract data for M and F groups
data_M <- subset(merged_RNAmDNA, Sex == "M")
data_F <- subset(merged_RNAmDNA, Sex == "F")

# Calculate correlation coefficients and p-values for M and F groups
cor_M <- cor.test(data_M$Ave_DNAm, data_M$ZNF727, method = "pearson")
cor_F <- cor.test(data_F$Ave_DNAm, data_F$ZNF727, method = "pearson")

# Create scatter plot for M and F groups
plot <- ggscatter(merged_RNAmDNA, x = "Ave_DNAm", y = "ZNF727", 
                  color = "Sex", palette = c("#E64B35B2", "#4DBBD5B2"),
                  add = "reg.line", conf.int = TRUE,
                  cor.coef = FALSE, p.value = FALSE) +
  xlab("Ave_DNAm") +
  ylab("Expression Level of ZNF727") +
  ggtitle("") +
  theme_bw()
print(plot)

plot_with_stats <- plot +
  geom_text(data = data_M, aes(label = paste("M: r =", round(cor_M$estimate, 2), "\n", "p =", format(cor_M$p.value, scientific = TRUE, digits = 2))),
            x = max(merged_RNAmDNA$Ave_DNAm), y = max(merged_RNAmDNA$ZNF727), hjust = 1, vjust = 3, color = "#4DBBD5B2") +
  geom_text(data = data_F, aes(label = paste("F: r =", round(cor_F$estimate, 2), "\n", "p =", format(cor_F$p.value, scientific = TRUE, digits = 2))),
            x = max(merged_RNAmDNA$Ave_DNAm), y = max(merged_RNAmDNA$ZNF727), hjust = 1, vjust = 1, color = "#E64B35B2") +
  theme(panel.grid = element_blank())

# Show the plot with correlation coefficients and p-values without grid lines
print(plot_with_stats)
ggsave("Cor_Ave_DNAm.pdf", plot_with_stats, width = 8, height = 6)



####——————————————————————————————————————————————————————————————————————————————————
####——————————————————————————————————————————————————————————————————————————————————
V4_sample = read.csv("ppmi_140_clinical_sPD_V04.csv", row.names = 1)
V6_sample = read.csv("ppmi_140_clinical_sPD_V06.csv", row.names = 1)
V8_sample = read.csv("ppmi_140_clinical_sPD_V08.csv", row.names = 1)
colnames(merged_RNAmDNA)
col_names <- colnames(merged_RNAmDNA)
new_col_names <- c("BL_Age", "BL_UPDRSpIA", "BL_UPDRSpIB", "BL_UPDRSpI", "BL_UPDRSpII", "BL_UPDRSpIII", "BL_MOCA", "BL_GDS")
colnames(merged_RNAmDNA)[12:19] <- new_col_names
print(colnames(merged_RNAmDNA))
library(dplyr)
matched_rows <- match(rownames(merged_RNAmDNA), rownames(V4_sample))
merged_RNAmDNAV4 <- bind_cols(merged_RNAmDNA, V4_sample[matched_rows, c("V4_UPDRSpIA", "V4_UPDRSpIB", "V4_UPDRSpI", "V4_UPDRSpII", "V4_UPDRSpIII", "V4_MOCA", "V4_GDS")])
print(merged_RNAmDNAV4)

matched_rows <- match(rownames(merged_RNAmDNAV4), rownames(V6_sample))
merged_RNAmDNAV6 <- bind_cols(merged_RNAmDNAV4, V6_sample[matched_rows, c("V6_UPDRSpIA", "V6_UPDRSpIB", "V6_UPDRSpI", "V6_UPDRSpII", "V6_UPDRSpIII", "V6_MOCA", "V6_GDS")])
print(merged_RNAmDNAV6)

matched_rows <- match(rownames(merged_RNAmDNAV6), rownames(V8_sample))
merged_RNAmDNAV8 <- bind_cols(merged_RNAmDNAV6, V8_sample[matched_rows, c("V8_UPDRSpIA", "V8_UPDRSpIB", "V8_UPDRSpI", "V8_UPDRSpII", "V8_UPDRSpIII", "V8_MOCA", "V8_GDS")])
print(merged_RNAmDNAV8)
merged_RNAmDNAclinical = merged_RNAmDNAV8
rm(merged_RNAmDNAV4)
rm(merged_RNAmDNAV6)
rm(merged_RNAmDNAV8)
colnames(merged_RNAmDNAclinical)
write.csv(merged_RNAmDNAclinical, "merged_RNAmDNAclinical.csv", row.names = TRUE)

missing_values_ratio <- colMeans(is.na(merged_RNAmDNAclinical))
print(missing_values_ratio)
column_indices <- 20:40
merged_RNAmDNAclinical <- merged_RNAmDNAclinical %>%
  mutate_at(column_indices, ~ifelse(is.na(.), mean(., na.rm = TRUE), .))
print(merged_RNAmDNAclinical)






###——————————————————————————————————————————————————————————————————————————————————
###——————————————————————————————————————————————————————————————————————————————————
install.packages("readxl")
library(readxl)
ImprintedGenes <- read_excel("Imprinted Genes_human.xlsx" )
ImprintedGenes = as.data.frame(ImprintedGenes)
class(ImprintedGenes)
ImprintedGenes$Gene

library(biomaRt)
library(dplyr)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_symbols <- unique(ImprintedGenes$Gene)
gene_info <- getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
                   filters = "external_gene_name",
                   values = gene_symbols,
                   mart = ensembl)

ImprintedGenesann <- merge(ImprintedGenes, gene_info, by.x = "Gene", by.y = "external_gene_name", all.x = TRUE)
ImprintedGenesann <- subset(ImprintedGenesann, !grepl("CHR", chromosome_name))
colnames(ImprintedGenesann) <- c("Gene", "Aliases", "Location", "Status", "Expressed Allele", "Chr", "Start", "End")
print(ImprintedGenesann)
write.csv(ImprintedGenesann, "ImprintedGenesann.csv", row.names = FALSE)

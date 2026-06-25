#2023-6-24

rm(list=ls())
options(stringsAsFactors = F) 


library(cit)
??cit

file_path <- "path/to/PDSex/eQTMs/DNAm_BL.csv"
DNAm_BL <- read.csv(file_path, row.names = 1)
targetCpG <- c("cg06098368", "cg03124146", "cg12479444", "cg20067334", "cg14285533", "cg26911611")
L <- DNAm_BL[rownames(DNAm_BL) %in% targetCpG, ]
L <- t(L)
row_names <- rownames(L)
new_row_names <- gsub("X", "", row_names)
new_row_names <- as.numeric(new_row_names)
rownames(L) <- new_row_names
class(L) #"matrix" "array" 

file_path <- "path/to/PDSex/eQTMs/eGene_BL.csv"
eGene_BL <- read.csv(file_path, row.names = 1)
G <- eGene_BL["ZNF727", , drop = FALSE] 
G <- t(G)
row_names <- rownames(G)
new_row_names <- gsub("X", "", row_names)
new_row_names <- as.numeric(new_row_names)
rownames(G) <- new_row_names
rows_match <- all(rownames(L) %in% rownames(G)) & all(rownames(G) %in% rownames(L)) #TRUE
class(G) #"matrix" "array" 

file_path <- "path/to/PPMI/Me_project 140/ppmi_140_sPD_BL/ppmi_140_clinical_sPD_BL.csv"
clinical_sPD_BL <- read.csv(file_path, row.names = 1) 
dim(clinical_sPD_BL)
rows <- rownames(clinical_sPD_BL) %in% rownames(G)
clinical_sPD_BL <- clinical_sPD_BL[rows, ]
rows_match <- identical(rownames(G), rownames(clinical_sPD_BL))
clinical_sPD_BL$Sex <- ifelse(clinical_sPD_BL$Sex == "F", 0, 1)
table(clinical_sPD_BL$Sex)
colnames(clinical_sPD_BL)
C <- clinical_sPD_BL[, "Age", drop = FALSE]
C <- as.matrix(C)


T_UPDRSpartIA_cp <- clinical_sPD_BL[, "MDS.UPDRS.part.IA", drop = FALSE]
T_UPDRSpartIB_cp <- clinical_sPD_BL[, "MDS.UPDRS.part.IB", drop = FALSE]
T_UPDRSpartI_cp <- clinical_sPD_BL[, "MDS.UPDRS.part.I", drop = FALSE]
T_UPDRSpartII_cp <- clinical_sPD_BL[, "MDS.UPDRS.part.II", drop = FALSE]
T_UPDRSpartIII_cp <- clinical_sPD_BL[, "MDS.UPDRS.part.III", drop = FALSE]
T_MOCA_cp <- clinical_sPD_BL[, "MOCA", drop = FALSE]
T_GDS_cp <- clinical_sPD_BL[, "GDS", drop = FALSE]
T_UPDRSpartIA_cp <- as.matrix(T_UPDRSpartIA_cp)
T_UPDRSpartIB_cp <- as.matrix(T_UPDRSpartIB_cp)
T_UPDRSpartI_cp <- as.matrix(T_UPDRSpartI_cp)
T_UPDRSpartII_cp <- as.matrix(T_UPDRSpartII_cp)
T_UPDRSpartIII_cp <- as.matrix(T_UPDRSpartIII_cp)
T_MOCA_cp <- as.matrix(T_MOCA_cp)
T_GDS_cp <- as.matrix(T_GDS_cp)

# Sample Size
ss = 199
n.perm = 5
perm.index = matrix(NA, nrow=ss, ncol=n.perm )
for( j in 1:ncol(perm.index) ) perm.index[, j] = sample( 1:ss )
cit_T_UPDRSpartIA_cp = cit.cp(L, G, T_UPDRSpartIA_cp, C)
cit_T_UPDRSpartIA_cp #0.9114708

cit_T_UPDRSpartIB_cp = cit.cp(L, G, T_UPDRSpartIB_cp, C)
cit_T_UPDRSpartIB_cp #0 0.8689022

cit_T_UPDRSpartI_cp = cit.cp(L, G, T_UPDRSpartI_cp, C)
cit_T_UPDRSpartI_cp #0 0.9348428

cit_T_UPDRSpartII_cp = cit.cp(L, G, T_UPDRSpartII_cp, C)
cit_T_UPDRSpartII_cp #0.4672128

cit_T_UPDRSpartIII_cp = cit.cp(L, G, T_UPDRSpartIII_cp, C)
cit_T_UPDRSpartIII_cp #0.6532614

cit_T_MOCA_cp = cit.cp(L, G, T_MOCA_cp, C)
cit_T_MOCA_cp #0.6921120

cit_T_GDS_cp = cit.cp(L, G, T_GDS_cp, C)
cit_T_GDS_cp #0.6184542


T_UPDRSpartIA_bp <- as.matrix(T_UPDRSpartIA_cp > median(T_UPDRSpartIA_cp)) + 0
T_UPDRSpartIB_bp <- as.matrix(T_UPDRSpartIB_cp > median(T_UPDRSpartIB_cp)) + 0
T_UPDRSpartI_bp <- as.matrix(T_UPDRSpartI_cp > median(T_UPDRSpartI_cp)) + 0
T_UPDRSpartII_bp <- as.matrix(T_UPDRSpartII_cp > median(T_UPDRSpartII_cp)) + 0
T_UPDRSpartIII_bp <- as.matrix(T_UPDRSpartIII_cp > median(T_UPDRSpartIII_cp)) + 0
T_MOCA_bp <- as.matrix(T_MOCA_cp > median(T_MOCA_cp)) + 0
T_GDS_bp <- as.matrix(T_GDS_cp > median(T_GDS_cp)) + 0

cit_T_UPDRSpartIA_bp = cit.bp(L, G, T_UPDRSpartIA_bp, C)
cit_T_UPDRSpartIA_bp #0.617824

cit_T_UPDRSpartIB_bp = cit.bp(L, G, T_UPDRSpartIB_bp, C)
cit_T_UPDRSpartIB_bp # 0.0948745

cit_T_UPDRSpartI_bp = cit.bp(L, G, T_UPDRSpartI_bp, C)
cit_T_UPDRSpartI_bp  # 0.4955826 

cit_T_UPDRSpartII_bp = cit.bp(L, G, T_UPDRSpartII_bp, C)
cit_T_UPDRSpartII_bp  #0.5568509

cit_T_UPDRSpartIII_bp = cit.bp(L, G, T_UPDRSpartIII_bp, C)
cit_T_UPDRSpartIII_bp #0.796

cit_T_MOCA_bp = cit.bp(L, G, T_MOCA_bp, C)
cit_T_MOCA_bp #0.7197691

cit_T_GDS_bp = cit.bp(L, G, T_GDS_bp, C)
cit_T_GDS_bp # 0.619

####__________________________________________________________________
file_path <- "path/to/PPMI/Me_project 140/ppmi_140_sPD_V08/ppmi_140_clinical_sPD_V08.csv"
clinical_sPD_V8 <- read.csv(file_path, row.names = 1) 
dim(clinical_sPD_V8)
rows <- rownames(clinical_sPD_V8) %in% rownames(G)
clinical_sPD_V8 <- clinical_sPD_V8[rows, ]
clinical_sPD_V8$Sex <- ifelse(clinical_sPD_V8$Sex == "F", 0, 1)
table(clinical_sPD_V8$Sex)
clinical_sPD_V8$MOCA <- as.numeric(as.character(clinical_sPD_V8$MOCA))
clinical_sPD_V8$MOCA[is.na(clinical_sPD_V8$MOCA)] <- mean(clinical_sPD_V8$MOCA, na.rm = TRUE)

colnames(clinical_sPD_V8)
T_V8_UPDRSpartIA_cp <- as.matrix(clinical_sPD_V8[, "MDS.UPDRS.part.IA", drop = FALSE])
T_V8_UPDRSpartIB_cp <- as.matrix(clinical_sPD_V8[, "MDS.UPDRS.part.IB", drop = FALSE])
T_V8_UPDRSpartI_cp <- as.matrix(clinical_sPD_V8[, "MDS.UPDRS.part.I", drop = FALSE])
T_V8_UPDRSpartII_cp <- as.matrix(clinical_sPD_V8[, "MDS.UPDRS.part.II", drop = FALSE])
T_V8_UPDRSpartIII_cp <- as.matrix(clinical_sPD_V8[, "MDS.UPDRS.part.III", drop = FALSE])
T_V8_MOCA_cp <- as.matrix(clinical_sPD_V8[, "MOCA", drop = FALSE])
T_V8_GDS_cp <- as.matrix(clinical_sPD_V8[, "GDS", drop = FALSE])
C_184 <- as.matrix(clinical_sPD_V8[, "Age", drop = FALSE])

common_rows <- rownames(L) %in% rownames(C_184)
L_184 <- L[common_rows, ]
G_184 <- G[common_rows, ]
G_184 = as.matrix(G_184)

cit_T_V8_UPDRSpartIA_cp = cit.cp(L_184, G_184, T_V8_UPDRSpartIA_cp, C_184)
cit_T_V8_UPDRSpartIA_cp #5.219311e-01

cit_T_V8_UPDRSpartIB_cp = cit.cp(L_184, G_184, T_V8_UPDRSpartIB_cp, C_184)
cit_T_V8_UPDRSpartIB_cp #9.972709e-01

cit_T_V8_UPDRSpartI_cp = cit.cp(L_184, G_184, T_V8_UPDRSpartI_cp, C_184)
cit_T_V8_UPDRSpartI_cp  #9.672702e-01

cit_T_V8_UPDRSpartII_cp = cit.cp(L_184, G_184, T_V8_UPDRSpartII_cp, C_184)
cit_T_V8_UPDRSpartII_cp  #5.044466e-01

cit_T_V8_UPDRSpartIII_cp = cit.cp(L_184, G_184, T_V8_UPDRSpartIII_cp, C_184)
cit_T_V8_UPDRSpartIII_cp #3.203118e-01

cit_T_V8_MOCA_cp = cit.cp(L_184, G_184, T_V8_MOCA_cp, C_184)
cit_T_V8_MOCA_cp #8.438409e-01

cit_T_V8_GDS_cp = cit.cp(L_184, G_184, T_V8_GDS_cp, C_184)
cit_T_V8_GDS_cp #8.509645e-01


########_______________________________________________________________________
T_V8_UPDRSpartIA_bp <- as.matrix(T_V8_UPDRSpartIA_cp > median(T_V8_UPDRSpartIA_cp)) + 0
T_V8_UPDRSpartIB_bp <- as.matrix(T_V8_UPDRSpartIB_cp > median(T_V8_UPDRSpartIB_cp)) + 0
T_V8_UPDRSpartI_bp <- as.matrix(T_V8_UPDRSpartI_cp > median(T_V8_UPDRSpartI_cp)) + 0
T_V8_UPDRSpartII_bp <- as.matrix(T_V8_UPDRSpartII_cp > median(T_V8_UPDRSpartII_cp)) + 0
T_V8_UPDRSpartIII_bp <- as.matrix(T_V8_UPDRSpartIII_cp > median(T_V8_UPDRSpartIII_cp)) + 0
T_V8_MOCA_bp <- as.matrix(T_V8_MOCA_cp > median(T_V8_MOCA_cp)) + 0
T_V8_GDS_bp <- as.matrix(T_V8_GDS_cp > median(T_V8_GDS_cp)) + 0
table(T_V8_UPDRSpartIA_bp)
table(T_V8_UPDRSpartIB_bp)
table(T_V8_UPDRSpartI_bp)
table(T_V8_UPDRSpartII_bp)
table(T_V8_UPDRSpartIII_bp)
table(T_V8_MOCA_bp)
table(T_V8_GDS_bp)

cit_T_V8_UPDRSpartIA_bp = cit.bp(L_184, G_184, T_V8_UPDRSpartIA_bp, C_184)
cit_T_V8_UPDRSpartIA_bp #0.5879036 

cit_T_V8_UPDRSpartIB_bp = cit.bp(L_184, G_184, T_V8_UPDRSpartIB_bp, C_184)
cit_T_V8_UPDRSpartIB_bp #0.4792122 

cit_T_V8_UPDRSpartI_bp = cit.bp(L_184, G_184, T_V8_UPDRSpartI_bp, C_184)
cit_T_V8_UPDRSpartI_bp  #0.141

cit_T_V8_UPDRSpartII_bp = cit.bp(L_184, G_184, T_V8_UPDRSpartII_bp, C_184)
cit_T_V8_UPDRSpartII_bp  #0.7880478 

cit_T_V8_UPDRSpartIII_bp = cit.bp(L_184, G_184, T_V8_UPDRSpartIII_bp, C_184)
cit_T_V8_UPDRSpartIII_bp #0.1444178

cit_T_V8_MOCA_bp = cit.bp(L_184, G_184, T_V8_MOCA_bp, C_184)
cit_T_V8_MOCA_bp #0.9737273

cit_T_V8_GDS_bp = cit.bp(L_184, G_184, T_V8_GDS_bp, C_184)
cit_T_V8_GDS_bp  #0.926

#————————————————————————————————————————————————————————————————————————————————
file_path <- "path/to/PPMI/Me_project 140/ppmi_140_sPD_V06/ppmi_140_clinical_sPD_V06.csv"
clinical_sPD_V6 <- read.csv(file_path, row.names = 1) 
dim(clinical_sPD_V6)
rows <- rownames(clinical_sPD_V6) %in% rownames(G)
clinical_sPD_V6 <- clinical_sPD_V6[rows, ]
clinical_sPD_V6$Sex <- ifelse(clinical_sPD_V6$Sex == "F", 0, 1)
table(clinical_sPD_V6$Sex)
clinical_sPD_V6$V6_MOCA <- as.numeric(as.character(clinical_sPD_V6$V6_MOCA))
clinical_sPD_V6$V6_MOCA[is.na(clinical_sPD_V6$V6_MOCA)] <- mean(clinical_sPD_V6$V6_MOCA, na.rm = TRUE)
clinical_sPD_V6$V6_GDS <- as.numeric(as.character(clinical_sPD_V6$V6_GDS))
clinical_sPD_V6$V6_GDS[is.na(clinical_sPD_V6$V6_GDS)] <- mean(clinical_sPD_V6$V6_GDS, na.rm = TRUE)

colnames(clinical_sPD_V6)
T_V6_UPDRSpartIA_cp <- as.matrix(clinical_sPD_V6[, "V6_UPDRSpIA", drop = FALSE])
T_V6_UPDRSpartIB_cp <- as.matrix(clinical_sPD_V6[, "V6_UPDRSpIB", drop = FALSE])
T_V6_UPDRSpartI_cp <- as.matrix(clinical_sPD_V6[, "V6_UPDRSpI", drop = FALSE])
T_V6_UPDRSpartII_cp <- as.matrix(clinical_sPD_V6[, "V6_UPDRSpII", drop = FALSE])
T_V6_UPDRSpartIII_cp <- as.matrix(clinical_sPD_V6[, "V6_UPDRSpIII", drop = FALSE])
T_V6_MOCA_cp <- as.matrix(clinical_sPD_V6[, "V6_MOCA", drop = FALSE])
T_V6_GDS_cp <- as.matrix(clinical_sPD_V6[, "V6_GDS", drop = FALSE])
C_185 <- as.matrix(clinical_sPD_V6[, "V6_Age", drop = FALSE])

common_rows <- rownames(L) %in% rownames(C_185)
L_185 <- L[common_rows, ]
G_185 <- G[common_rows, ]
G_185 = as.matrix(G_185)

library(cit)
cit_T_V6_UPDRSpartIA_cp = cit.cp(L_185, G_185, T_V6_UPDRSpartIA_cp, C_185)
cit_T_V6_UPDRSpartIA_cp #9.659812e-01

cit_T_V6_UPDRSpartIB_cp = cit.cp(L_185, G_185, T_V6_UPDRSpartIB_cp, C_185)
cit_T_V6_UPDRSpartIB_cp #8.742664e-01

cit_T_V6_UPDRSpartI_cp = cit.cp(L_185, G_185, T_V6_UPDRSpartI_cp, C_185)
cit_T_V6_UPDRSpartI_cp  #9.506307e-01

cit_T_V6_UPDRSpartII_cp = cit.cp(L_185, G_185, T_V6_UPDRSpartII_cp, C_185)
cit_T_V6_UPDRSpartII_cp  #3.148012e-01

cit_T_V6_UPDRSpartIII_cp = cit.cp(L_185, G_185, T_V6_UPDRSpartIII_cp, C_185)
cit_T_V6_UPDRSpartIII_cp #2.447394e-01

cit_T_V6_MOCA_cp = cit.cp(L_185, G_185, T_V6_MOCA_cp, C_185)
cit_T_V6_MOCA_cp #9.270694e-01

cit_T_V6_GDS_cp = cit.cp(L_185, G_185, T_V6_GDS_cp, C_185)
cit_T_V6_GDS_cp #9.953114e-01


T_V6_UPDRSpartIA_bp <- as.matrix(T_V6_UPDRSpartIA_cp > median(T_V6_UPDRSpartIA_cp)) + 0
T_V6_UPDRSpartIB_bp <- as.matrix(T_V6_UPDRSpartIB_cp > median(T_V6_UPDRSpartIB_cp)) + 0
T_V6_UPDRSpartI_bp <- as.matrix(T_V6_UPDRSpartI_cp > median(T_V6_UPDRSpartI_cp)) + 0
T_V6_UPDRSpartII_bp <- as.matrix(T_V6_UPDRSpartII_cp > median(T_V6_UPDRSpartII_cp)) + 0
T_V6_UPDRSpartIII_bp <- as.matrix(T_V6_UPDRSpartIII_cp > median(T_V6_UPDRSpartIII_cp)) + 0
T_V6_MOCA_bp <- as.matrix(T_V6_MOCA_cp > median(T_V6_MOCA_cp)) + 0
T_V6_GDS_bp <- as.matrix(T_V6_GDS_cp > median(T_V6_GDS_cp)) + 0
table(T_V6_UPDRSpartIA_bp)
table(T_V6_UPDRSpartIB_bp)
table(T_V6_UPDRSpartI_bp)
table(T_V6_UPDRSpartII_bp)
table(T_V6_UPDRSpartIII_bp)
table(T_V6_MOCA_bp)
table(T_V6_GDS_bp)

cit_T_V6_UPDRSpartIA_bp = cit.bp(L_185, G_185, T_V6_UPDRSpartIA_bp, C_185)
cit_T_V6_UPDRSpartIA_bp #0.9865293

cit_T_V6_UPDRSpartIB_bp = cit.bp(L_185, G_185, T_V6_UPDRSpartIB_bp, C_185)
cit_T_V6_UPDRSpartIB_bp #0.7150677 

cit_T_V6_UPDRSpartI_bp = cit.bp(L_185, G_185, T_V6_UPDRSpartI_bp, C_185)
cit_T_V6_UPDRSpartI_bp  #0.9716311 

cit_T_V6_UPDRSpartII_bp = cit.bp(L_185, G_185, T_V6_UPDRSpartII_bp, C_185)
cit_T_V6_UPDRSpartII_bp  #0.5983641  

cit_T_V6_UPDRSpartIII_bp = cit.bp(L_185, G_185, T_V6_UPDRSpartIII_bp, C_185)
cit_T_V6_UPDRSpartIII_bp #0.3364889 

cit_T_V6_MOCA_bp = cit.bp(L_185, G_185, T_V6_MOCA_bp, C_185)
cit_T_V6_MOCA_bp #0.6984354 

cit_T_V6_GDS_bp = cit.bp(L_185, G_185, T_V6_GDS_bp, C_185)
cit_T_V6_GDS_bp #0.943 

#——————————————————————————————————————————————————————————————————————————————————————



#————————————————————————————————————————————————————————————————————————————————
file_path <- "path/to/PPMI/Me_project 140/ppmi_140_sPD_V04/ppmi_140_clinical_sPD_V04.csv"
clinical_sPD_V4 <- read.csv(file_path, row.names = 1) 
dim(clinical_sPD_V4)
rows <- rownames(clinical_sPD_V4) %in% rownames(G)
clinical_sPD_V4 <- clinical_sPD_V4[rows, ]
clinical_sPD_V4$Sex <- ifelse(clinical_sPD_V4$Sex == "F", 0, 1)
table(clinical_sPD_V4$Sex)
clinical_sPD_V4$V4_MOCA <- as.numeric(as.character(clinical_sPD_V4$V4_MOCA))
clinical_sPD_V4$V4_MOCA[is.na(clinical_sPD_V4$V4_MOCA)] <- mean(clinical_sPD_V4$V4_MOCA, na.rm = TRUE)
clinical_sPD_V4$V4_GDS <- as.numeric(as.character(clinical_sPD_V4$V4_GDS))
clinical_sPD_V4$V4_GDS[is.na(clinical_sPD_V4$V4_GDS)] <- mean(clinical_sPD_V4$V4_GDS, na.rm = TRUE)

colnames(clinical_sPD_V4)
T_V4_UPDRSpartIA_cp <- as.matrix(clinical_sPD_V4[, "V4_UPDRSpIA", drop = FALSE])
T_V4_UPDRSpartIB_cp <- as.matrix(clinical_sPD_V4[, "V4_UPDRSpIB", drop = FALSE])
T_V4_UPDRSpartI_cp <- as.matrix(clinical_sPD_V4[, "V4_UPDRSpI", drop = FALSE])
T_V4_UPDRSpartII_cp <- as.matrix(clinical_sPD_V4[, "V4_UPDRSpII", drop = FALSE])
T_V4_UPDRSpartIII_cp <- as.matrix(clinical_sPD_V4[, "V4_UPDRSpIII", drop = FALSE])
T_V4_MOCA_cp <- as.matrix(clinical_sPD_V4[, "V4_MOCA", drop = FALSE])
T_V4_GDS_cp <- as.matrix(clinical_sPD_V4[, "V4_GDS", drop = FALSE])
C_177 <- as.matrix(clinical_sPD_V4[, "V4_Age", drop = FALSE])

common_rows <- rownames(L) %in% rownames(C_177)
L_177 <- L[common_rows, ]
G_177 <- G[common_rows, ]
G_177 = as.matrix(G_177)

library(cit)
cit_T_V4_UPDRSpartIA_cp = cit.cp(L_177, G_177, T_V4_UPDRSpartIA_cp, C_177)
cit_T_V4_UPDRSpartIA_cp #9.159493e-01

cit_T_V4_UPDRSpartIB_cp = cit.cp(L_177, G_177, T_V4_UPDRSpartIB_cp, C_177)
cit_T_V4_UPDRSpartIB_cp #6.114946e-01

cit_T_V4_UPDRSpartI_cp = cit.cp(L_177, G_177, T_V4_UPDRSpartI_cp, C_177)
cit_T_V4_UPDRSpartI_cp  #7.340702e-01 

cit_T_V4_UPDRSpartII_cp = cit.cp(L_177, G_177, T_V4_UPDRSpartII_cp, C_177)
cit_T_V4_UPDRSpartII_cp  #4.034926e-01

cit_T_V4_UPDRSpartIII_cp = cit.cp(L_177, G_177, T_V4_UPDRSpartIII_cp, C_177)
cit_T_V4_UPDRSpartIII_cp #6.341087e-01

cit_T_V4_MOCA_cp = cit.cp(L_177, G_177, T_V4_MOCA_cp, C_177)
cit_T_V4_MOCA_cp #4.140319e-01  

cit_T_V4_GDS_cp = cit.cp(L_177, G_177, T_V4_GDS_cp, C_177)
cit_T_V4_GDS_cp # 5.921412e-01


T_V4_UPDRSpartIA_bp <- as.matrix(T_V4_UPDRSpartIA_cp > median(T_V4_UPDRSpartIA_cp)) + 0
T_V4_UPDRSpartIB_bp <- as.matrix(T_V4_UPDRSpartIB_cp > median(T_V4_UPDRSpartIB_cp)) + 0
T_V4_UPDRSpartI_bp <- as.matrix(T_V4_UPDRSpartI_cp > median(T_V4_UPDRSpartI_cp)) + 0
T_V4_UPDRSpartII_bp <- as.matrix(T_V4_UPDRSpartII_cp > median(T_V4_UPDRSpartII_cp)) + 0
T_V4_UPDRSpartIII_bp <- as.matrix(T_V4_UPDRSpartIII_cp > median(T_V4_UPDRSpartIII_cp)) + 0
T_V4_MOCA_bp <- as.matrix(T_V4_MOCA_cp > median(T_V4_MOCA_cp)) + 0
T_V4_GDS_bp <- as.matrix(T_V4_GDS_cp > median(T_V4_GDS_cp)) + 0
table(T_V4_UPDRSpartIA_bp)
table(T_V4_UPDRSpartIB_bp)
table(T_V4_UPDRSpartI_bp)
table(T_V4_UPDRSpartII_bp)
table(T_V4_UPDRSpartIII_bp)
table(T_V4_MOCA_bp)
table(T_V4_GDS_bp)

cit_T_V4_UPDRSpartIA_bp = cit.bp(L_177, G_177, T_V4_UPDRSpartIA_bp, C_177)
cit_T_V4_UPDRSpartIA_bp #0.695

cit_T_V4_UPDRSpartIB_bp = cit.bp(L_177, G_177, T_V4_UPDRSpartIB_bp, C_177)
cit_T_V4_UPDRSpartIB_bp #0.6565695

cit_T_V4_UPDRSpartI_bp = cit.bp(L_177, G_177, T_V4_UPDRSpartI_bp, C_177)
cit_T_V4_UPDRSpartI_bp  #0.7549683 

cit_T_V4_UPDRSpartII_bp = cit.bp(L_177, G_177, T_V4_UPDRSpartII_bp, C_177)
cit_T_V4_UPDRSpartII_bp  #0.8440278 

cit_T_V4_UPDRSpartIII_bp = cit.bp(L_177, G_177, T_V4_UPDRSpartIII_bp, C_177)
cit_T_V4_UPDRSpartIII_bp #0.7389232 

cit_T_V4_MOCA_bp = cit.bp(L_177, G_177, T_V4_MOCA_bp, C_177)
cit_T_V4_MOCA_bp #0.7942534 

cit_T_V4_GDS_bp = cit.bp(L_177, G_177, T_V4_GDS_bp, C_177)
cit_T_V4_GDS_bp #0.649  

#——————————————————————————————————————————————————————————————————————————————————————



L_all = t(DNAm_BL)
class(L_all)
rownames(L_all) <- gsub("X", "", rownames(L_all))
common_rows <- intersect(rownames(L_all), rownames(clinical_sPD_V8))
L_all <- L_all[common_rows, ]


citall_T_V8_UPDRSpartIA_bp = cit.bp(L_all, G_184, T_V8_UPDRSpartIA_bp, C_184)
citall_T_V8_UPDRSpartIA_bp #0.5879036

citall_T_V8_UPDRSpartIB_bp = cit.bp(L_all, G_184, T_V8_UPDRSpartIB_bp, C_184)
citall_T_V8_UPDRSpartIB_bp #0.4792122

citall_T_UPDRSpartI_bp = cit.bp(L_all, G_184, T_V8_UPDRSpartI_bp, C_184)
citall_T_UPDRSpartI_bp  #0.145

citall_T_V8_UPDRSpartII_bp = cit.bp(L_all, G_184, T_V8_UPDRSpartII_bp, C_184)
citall_T_V8_UPDRSpartII_bp  #0.7880478

citall_T_V8_UPDRSpartIII_bp = cit.bp(L_all, G_184, T_V8_UPDRSpartIII_bp, C_184)
citall_T_V8_UPDRSpartIII_bp #0.1444178

citall_T_V8_MOCA_bp = cit.bp(L_all, G_184, T_V8_MOCA_bp, C_184)
citall_T_V8_MOCA_bp #0.9737273

citall_T_V8_GDS_bp = cit.bp(L_all, G_184, T_V8_GDS_bp, C_184)
citall_T_V8_GDS_bp  #0.924

#——————————————————————————————————————————————————————————————————————
#BL-CP
cit_T_UPDRSpartIA_cp #0.9114708
cit_T_UPDRSpartIB_cp #0 0.8689022
cit_T_UPDRSpartI_cp #0 0.9348428
cit_T_UPDRSpartII_cp #0.4672128
cit_T_UPDRSpartIII_cp #0.6532614
cit_T_MOCA_cp #0.6921120
cit_T_GDS_cp #0.6184542
cit_T_UPDRSpartIA_cp <- as.data.frame(t(cit_T_UPDRSpartIA_cp))
cit_T_UPDRSpartIB_cp <- as.data.frame(t(cit_T_UPDRSpartIB_cp))
cit_T_UPDRSpartI_cp <- as.data.frame(t(cit_T_UPDRSpartI_cp))
cit_T_UPDRSpartII_cp <- as.data.frame(t(cit_T_UPDRSpartII_cp))
cit_T_UPDRSpartIII_cp <- as.data.frame(t(cit_T_UPDRSpartIII_cp))
cit_T_MOCA_cp <- as.data.frame(t(cit_T_MOCA_cp))
cit_T_GDS_cp <- as.data.frame(t(cit_T_GDS_cp))
write.csv(cit_T_UPDRSpartIA_cp, "cit_T_UPDRSpartIA_cp.csv", row.names = FALSE)
write.csv(cit_T_UPDRSpartIB_cp, "cit_T_UPDRSpartIB_cp.csv", row.names = FALSE)
write.csv(cit_T_UPDRSpartI_cp, "cit_T_UPDRSpartI_cp.csv", row.names = FALSE)
write.csv(cit_T_UPDRSpartII_cp, "cit_T_UPDRSpartII_cp.csv", row.names = FALSE)
write.csv(cit_T_UPDRSpartIII_cp, "cit_T_UPDRSpartIII_cp.csv", row.names = FALSE)
write.csv(cit_T_MOCA_cp, "cit_T_MOCA_cp.csv", row.names = FALSE)
write.csv(cit_T_GDS_cp, "cit_T_GDS_cp.csv", row.names = FALSE)

#BL-BP
cit_T_UPDRSpartIA_bp #0.617824
cit_T_UPDRSpartIB_bp #0.0948745
cit_T_UPDRSpartI_bp  #0.4955826 
cit_T_UPDRSpartII_bp  #0.5568509 
cit_T_UPDRSpartIII_bp #0.796 
cit_T_MOCA_bp #0.7197691 
cit_T_GDS_bp  #0.619
class(cit_T_UPDRSpartIA_bp)

write.csv(cit_T_UPDRSpartIA_bp, "cit_T_UPDRSpartIA_bp.csv", row.names = FALSE)
write.csv(cit_T_UPDRSpartIB_bp, "cit_T_UPDRSpartIB_bp.csv", row.names = FALSE)
write.csv(cit_T_UPDRSpartI_bp, "cit_T_UPDRSpartI_bp.csv", row.names = FALSE)
write.csv(cit_T_UPDRSpartII_bp, "cit_T_UPDRSpartII_bp.csv", row.names = FALSE)
write.csv(cit_T_UPDRSpartIII_bp, "cit_T_UPDRSpartIII_bp.csv", row.names = FALSE)
write.csv(cit_T_MOCA_bp, "cit_T_MOCA_bp.csv", row.names = FALSE)
write.csv(cit_T_GDS_bp, "cit_T_GDS_bp.csv", row.names = FALSE)

#________________________________________________________________
#V8-CP
cit_T_V8_UPDRSpartIA_cp #5.219311e-01
cit_T_V8_UPDRSpartIB_cp #9.972709e-01
cit_T_V8_UPDRSpartI_cp  #9.672702e-01 
cit_T_V8_UPDRSpartII_cp  #5.044466e-01
cit_T_V8_UPDRSpartIII_cp #3.203118e-01
cit_T_V8_MOCA_cp #8.438409e-01
cit_T_V8_GDS_cp #8.509645e-01
cit_T_V8_UPDRSpartIA_cp <- as.data.frame(t(cit_T_V8_UPDRSpartIA_cp))
cit_T_V8_UPDRSpartIB_cp <- as.data.frame(t(cit_T_V8_UPDRSpartIB_cp))
cit_T_V8_UPDRSpartI_cp <- as.data.frame(t(cit_T_V8_UPDRSpartI_cp))
cit_T_V8_UPDRSpartII_cp <- as.data.frame(t(cit_T_V8_UPDRSpartII_cp))
cit_T_V8_UPDRSpartIII_cp <- as.data.frame(t(cit_T_V8_UPDRSpartIII_cp))
cit_T_V8_MOCA_cp <- as.data.frame(t(cit_T_V8_MOCA_cp))
cit_T_V8_GDS_cp <- as.data.frame(t(cit_T_V8_GDS_cp))
write.csv(cit_T_V8_UPDRSpartIA_cp, "cit_T_V8_UPDRSpartIA_cp.csv", row.names = FALSE)
write.csv(cit_T_V8_UPDRSpartIB_cp, "cit_T_V8_UPDRSpartIB_cp.csv", row.names = FALSE)
write.csv(cit_T_V8_UPDRSpartI_cp, "cit_T_V8_UPDRSpartI_cp.csv", row.names = FALSE)
write.csv(cit_T_V8_UPDRSpartII_cp, "cit_T_V8_UPDRSpartII_cp.csv", row.names = FALSE)
write.csv(cit_T_V8_UPDRSpartIII_cp, "cit_T_V8_UPDRSpartIII_cp.csv", row.names = FALSE)
write.csv(cit_T_V8_MOCA_cp, "cit_T_V8_MOCA_cp.csv", row.names = FALSE)
write.csv(cit_T_V8_GDS_cp, "cit_T_V8_GDS_cp.csv", row.names = FALSE)


#V8-BP
cit_T_V8_UPDRSpartIA_bp #0.5879036 
cit_T_V8_UPDRSpartIB_bp #0.4792122 
cit_T_V8_UPDRSpartI_bp  #0.141
cit_T_V8_UPDRSpartII_bp  #0.7880478 
cit_T_V8_UPDRSpartIII_bp #0.1444178
cit_T_V8_MOCA_bp #0.9737273
cit_T_V8_GDS_bp  #0.926
class(cit_T_V8_UPDRSpartIA_bp)

write.csv(cit_T_V8_UPDRSpartIA_bp, "cit_T_V8_UPDRSpartIA_bp.csv", row.names = FALSE)
write.csv(cit_T_V8_UPDRSpartIB_bp, "cit_T_V8_UPDRSpartIB_bp.csv", row.names = FALSE)
write.csv(cit_T_V8_UPDRSpartI_bp, "cit_T_V8_UPDRSpartI_bp.csv", row.names = FALSE)
write.csv(cit_T_V8_UPDRSpartII_bp, "cit_T_V8_UPDRSpartII_bp.csv", row.names = FALSE)
write.csv(cit_T_V8_UPDRSpartIII_bp, "cit_T_V8_UPDRSpartIII_bp.csv", row.names = FALSE)
write.csv(cit_T_V8_MOCA_bp, "cit_T_V8_MOCA_bp.csv", row.names = FALSE)
write.csv(cit_T_V8_GDS_bp, "cit_T_V8_GDS_bp.csv", row.names = FALSE)



#________________________________________________________________
#V6-CP
cit_T_V6_UPDRSpartIA_cp <- as.data.frame(t(cit_T_V6_UPDRSpartIA_cp))
cit_T_V6_UPDRSpartIB_cp <- as.data.frame(t(cit_T_V6_UPDRSpartIB_cp))
cit_T_V6_UPDRSpartI_cp <- as.data.frame(t(cit_T_V6_UPDRSpartI_cp))
cit_T_V6_UPDRSpartII_cp <- as.data.frame(t(cit_T_V6_UPDRSpartII_cp))
cit_T_V6_UPDRSpartIII_cp <- as.data.frame(t(cit_T_V6_UPDRSpartIII_cp))
cit_T_V6_MOCA_cp <- as.data.frame(t(cit_T_V6_MOCA_cp))
cit_T_V6_GDS_cp <- as.data.frame(t(cit_T_V6_GDS_cp))
write.csv(cit_T_V6_UPDRSpartIA_cp, "cit_T_V6_UPDRSpartIA_cp.csv", row.names = FALSE)
write.csv(cit_T_V6_UPDRSpartIB_cp, "cit_T_V6_UPDRSpartIB_cp.csv", row.names = FALSE)
write.csv(cit_T_V6_UPDRSpartI_cp, "cit_T_V6_UPDRSpartI_cp.csv", row.names = FALSE)
write.csv(cit_T_V6_UPDRSpartII_cp, "cit_T_V6_UPDRSpartII_cp.csv", row.names = FALSE)
write.csv(cit_T_V6_UPDRSpartIII_cp, "cit_T_V6_UPDRSpartIII_cp.csv", row.names = FALSE)
write.csv(cit_T_V6_MOCA_cp, "cit_T_V6_MOCA_cp.csv", row.names = FALSE)
write.csv(cit_T_V6_GDS_cp, "cit_T_V6_GDS_cp.csv", row.names = FALSE)


#V6-bp
write.csv(cit_T_V6_UPDRSpartIA_bp, "cit_T_V6_UPDRSpartIA_bp.csv", row.names = FALSE)
write.csv(cit_T_V6_UPDRSpartIB_bp, "cit_T_V6_UPDRSpartIB_bp.csv", row.names = FALSE)
write.csv(cit_T_V6_UPDRSpartI_bp, "cit_T_V6_UPDRSpartI_bp.csv", row.names = FALSE)
write.csv(cit_T_V6_UPDRSpartII_bp, "cit_T_V6_UPDRSpartII_bp.csv", row.names = FALSE)
write.csv(cit_T_V6_UPDRSpartIII_bp, "cit_T_V6_UPDRSpartIII_bp.csv", row.names = FALSE)
write.csv(cit_T_V6_MOCA_bp, "cit_T_V6_MOCA_bp.csv", row.names = FALSE)
write.csv(cit_T_V6_GDS_bp, "cit_T_V6_GDS_bp.csv", row.names = FALSE)


#________________________________________________________________
#V4-CP
cit_T_V4_UPDRSpartIA_cp <- as.data.frame(t(cit_T_V4_UPDRSpartIA_cp))
cit_T_V4_UPDRSpartIB_cp <- as.data.frame(t(cit_T_V4_UPDRSpartIB_cp))
cit_T_V4_UPDRSpartI_cp <- as.data.frame(t(cit_T_V4_UPDRSpartI_cp))
cit_T_V4_UPDRSpartII_cp <- as.data.frame(t(cit_T_V4_UPDRSpartII_cp))
cit_T_V4_UPDRSpartIII_cp <- as.data.frame(t(cit_T_V4_UPDRSpartIII_cp))
cit_T_V4_MOCA_cp <- as.data.frame(t(cit_T_V4_MOCA_cp))
cit_T_V4_GDS_cp <- as.data.frame(t(cit_T_V4_GDS_cp))
write.csv(cit_T_V4_UPDRSpartIA_cp, "cit_T_V4_UPDRSpartIA_cp.csv", row.names = FALSE)
write.csv(cit_T_V4_UPDRSpartIB_cp, "cit_T_V4_UPDRSpartIB_cp.csv", row.names = FALSE)
write.csv(cit_T_V4_UPDRSpartI_cp, "cit_T_V4_UPDRSpartI_cp.csv", row.names = FALSE)
write.csv(cit_T_V4_UPDRSpartII_cp, "cit_T_V4_UPDRSpartII_cp.csv", row.names = FALSE)
write.csv(cit_T_V4_UPDRSpartIII_cp, "cit_T_V4_UPDRSpartIII_cp.csv", row.names = FALSE)
write.csv(cit_T_V4_MOCA_cp, "cit_T_V4_MOCA_cp.csv", row.names = FALSE)
write.csv(cit_T_V4_GDS_cp, "cit_T_V4_GDS_cp.csv", row.names = FALSE)


#V4-bp
write.csv(cit_T_V4_UPDRSpartIA_bp, "cit_T_V4_UPDRSpartIA_bp.csv", row.names = FALSE)
write.csv(cit_T_V4_UPDRSpartIB_bp, "cit_T_V4_UPDRSpartIB_bp.csv", row.names = FALSE)
write.csv(cit_T_V4_UPDRSpartI_bp, "cit_T_V4_UPDRSpartI_bp.csv", row.names = FALSE)
write.csv(cit_T_V4_UPDRSpartII_bp, "cit_T_V4_UPDRSpartII_bp.csv", row.names = FALSE)
write.csv(cit_T_V4_UPDRSpartIII_bp, "cit_T_V4_UPDRSpartIII_bp.csv", row.names = FALSE)
write.csv(cit_T_V4_MOCA_bp, "cit_T_V4_MOCA_bp.csv", row.names = FALSE)
write.csv(cit_T_V4_GDS_bp, "cit_T_V4_GDS_bp.csv", row.names = FALSE)


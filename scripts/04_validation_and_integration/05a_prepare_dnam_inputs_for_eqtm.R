
##2023-6-20
rm(list=ls())
options(stringsAsFactors = F) 
load(file = 'ppmi140_champ_myRefbase_BL_PD213.Rdata')
beta = myRefbasePD$CorrectedBeta

file_path <- "path/to/PDSex/eQTMs/BL_1985DMP.csv"
DMP_BL <- read.csv(file_path)

class(beta) #"matrix" "array" 
DNAm <- beta[rownames(beta) %in% DMP_BL$X, ] 
colnames(DNAm) <- pD$PATNO

file_path <- "path/to/PPMI/RNA/PPMI_RNAseq_IR3_Analysis/R featurecount BL/metaDataBL.csv"
RNAsample <- read.csv(file_path)
DNAmsample = pD
cmsample <- intersect(RNAsample$PATNO, DNAmsample$PATNO)
DNAm_BL <- DNAm[, as.character(cmsample)]
DNAm_BL1 <- DNAm_BL[, order(colnames(DNAm_BL))]
DNAm_BL = DNAm_BL1
rm(DNAm_BL1)
write.csv(DNAm_BL, "DNAm_BL.csv", row.names = TRUE)

cpgloc <- DMP_BL[, c("X", "M_to_F.CHR", "M_to_F.MAPINFO")]
colnames(cpgloc) <- c("cpgid", "chr", "pos")
cpgloc <- cpgloc[match(rownames(DNAm_BL), cpgloc$cpgid), ]
write.csv(cpgloc, "cpgloc_BL.csv", row.names = F) 

cova_BL <- read.csv("ppmi_140_clinical_sPD_BL.csv", header = TRUE)
cova_BL <- cova_BL[, c("PATNO", "Age")]
cova_BL <- t(cova_BL)
colnames(cova_BL) <- as.character(cova_BL[1, ])
cova_BL <- cova_BL[, as.character(cmsample)]
cova_BL <- cova_BL[, order(colnames(cova_BL))]
write.csv(cova_BL, "cova_BL.csv", row.names = T) 







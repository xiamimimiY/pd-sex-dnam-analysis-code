cmPD = read.csv("cmPD.csv")
rownames(cmPD)=cmPD[,1]
uniquePD = read.csv("uniquePD.csv")
rownames(uniquePD)=uniquePD[,1]

library(ChAMP)

PDDMP <- read.csv("myDMPPD331_0.05_bacon.csv")
rownames(PDDMP)=PDDMP[,1]
colnames(PDDMP)[1] = 'CpGs'
colnames(cmPD)[1] = 'CpGs'
colnames(uniquePD)[1] = 'CpGs'
cmPDann = merge(cmPD,PDDMP,by="CpGs")
uniquePDann = merge(uniquePD,PDDMP,by="CpGs")
write.table(cmPDann, file = "cmPDann.csv",sep = ",", row.names=FALSE,col.names=TRUE)
write.table(uniquePDann, file = "uniquePDann.csv",sep = ",", row.names=FALSE,col.names=TRUE)

######
colnames(uniquePDann)
selected_cols <- c("CpGs", "Effect", "StdErr", "P.value", "Direction", "F_to_M.CHR", "F_to_M.MAPINFO", "F_to_M.gene", 
                   "F_to_M.feature", "F_to_M.cgi")

uniquePDann0.05 <- uniquePDann[selected_cols]
uniquePDann0.05$Classification <- "PD-unique-DMPs"
cmPDann0.05  <- cmPDann[selected_cols]
cmPDann0.05$Classification <- "PD-HC-cmDMPs"
PDmetaAnn <- rbind(uniquePDann0.05, cmPDann0.05)
colnames(PDmetaAnn)
colnames(PDmetaAnn) <- sub("^F_to_M\\.", "", colnames(PDmetaAnn))
colnames(PDmetaAnn) <- gsub("MAPINFO", "Position", colnames(PDmetaAnn))
PDmetaAnn <- PDmetaAnn[order(PDmetaAnn$P.value), ]
table(PDmetaAnn$Classification)  #PD-HC-cmDMPs PD-unique-DMPs  826           2199 
write.csv(PDmetaAnn, file = "0.05PDmetaAnn.csv", row.names = FALSE)

HCmeta <- read.csv("0.05metaHC.csv")
file_path <- "path/to/PDSex/blood_111629_571/ALL571/myDMPHC236_0.05_bacon.csv"
HCAnn <- read.csv(file_path)
colnames(HCAnn)[1] <- "CpGs"
colnames(HCmeta)[1] <- "CpGs"
matched_rows <- HCAnn[HCAnn$CpGs %in% HCmeta$CpGs, ]
colnames(matched_rows)
selected_cols <- c("CpGs", "F_to_M.CHR", "F_to_M.MAPINFO", "F_to_M.Strand", "F_to_M.Type", "F_to_M.gene", "F_to_M.feature", "F_to_M.cgi")
HCmetaAnn <- merge(HCmeta, matched_rows[selected_cols], by = "CpGs")
HCmetaAnn <- HCmetaAnn[, -c(2, 3)]
colnames(HCmetaAnn) <- sub("^F_to_M\\.", "", colnames(HCmetaAnn))
HCmetaAnn <- HCmetaAnn[, -c(8, 9)]
colnames(HCmetaAnn) <- gsub("MAPINFO", "Position", colnames(HCmetaAnn))
HCmetaAnn <- HCmetaAnn[order(HCmetaAnn$P.value), ] 
write.csv(HCmetaAnn, file = "0.05HCmetaAnn.csv", row.names = FALSE)

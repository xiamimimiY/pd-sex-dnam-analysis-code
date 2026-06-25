cmPD0.05 <- read.csv("cmPD0.05.csv") #825
uniquePD0.05 <- read.csv("uniquePD0.05.csv") #2199
PPMI_BL_0.05 <- read.csv("myDMP0.05_PD310_bacon.csv",header = TRUE, row.names = 1) #12745
PPMI_V4_0.05 <- read.csv("myDMP0.05_V4_PD190_bacon.csv") #186657
PPMI_V6_0.05 <- read.csv("myDMP0.05_V6_PD195_bacon.csv") #262940
PPMI_V8_0.05 <- read.csv("myDMP0.05_V8_PD197_bacon.csv") #309014
saliva0.05 <- read.csv("myDMPPD128_0.05_bacon.csv") #2674
PPMI_BL_0.05 <- subset(PPMI_BL_0.05, fdr.bacon < 0.05) #2402
PPMI_V4_0.05 <- subset(PPMI_V4_0.05, fdr.bacon < 0.05) #430
PPMI_V6_0.05 <- subset(PPMI_V6_0.05, fdr.bacon < 0.05) #291
PPMI_V8_0.05 <- subset(PPMI_V8_0.05, fdr.bacon < 0.05) #384
saliva0.05 <- subset(saliva0.05, fdr.bacon < 0.05) #2674
PPMI_BL_0.05$MarkerName <- row.names(PPMI_BL_0.05)
PPMI_V4_0.05$MarkerName <- row.names(PPMI_V4_0.05)
PPMI_V6_0.05$MarkerName <- row.names(PPMI_V6_0.05)
PPMI_V8_0.05$MarkerName <- row.names(PPMI_V8_0.05)
saliva0.05$MarkerName <- row.names(saliva0.05)

PPMI_common <- Reduce(intersect, list(PPMI_BL_0.05$MarkerName,PPMI_V4_0.05$MarkerName, PPMI_V6_0.05$MarkerName,PPMI_V8_0.05$MarkerName)) #119

meta_PPMI_BL<- intersect(uniquePD0.05$MarkerName,PPMI_BL_0.05$MarkerName) #536
meta_PPMI_V4<- intersect(uniquePD0.05$MarkerName,PPMI_V4_0.05$MarkerName) #53
meta_PPMI_V6<- intersect(uniquePD0.05$MarkerName,PPMI_V6_0.05$MarkerName) #24
meta_PPMI_V8<- intersect(uniquePD0.05$MarkerName,PPMI_V8_0.05$MarkerName) #53
meta_PPMIall <- intersect(uniquePD0.05$MarkerName,PPMI_common) #7
meta_saliva <- intersect(uniquePD0.05$MarkerName,saliva0.05$MarkerName) #536
allcommon <- intersect(metaU_PPMI,metaU_saliva) #4
rm(meta_PPMIall_ann)
AnnBL_meta_PPMIall = PPMI_BL_0.05[PPMI_BL_0.05$MarkerName %in% meta_PPMIall,]
AnnV4_meta_PPMIall = PPMI_V4_0.05[PPMI_V4_0.05$MarkerName %in% meta_PPMIall,]
AnnV6_meta_PPMIall = PPMI_V6_0.05[PPMI_V6_0.05$MarkerName %in% meta_PPMIall,]
AnnV8_meta_PPMIall = PPMI_V8_0.05[PPMI_V8_0.05$MarkerName %in% meta_PPMIall,]
Annsaliva_meta_PPMIall = saliva0.05[saliva0.05$MarkerName %in% allcommon,]
write.table(AnnBL_meta_PPMIall, file = "AnnBL_meta_PPMIall.csv", sep = ",", row.names=F,col.names=TRUE)
write.table(AnnV4_meta_PPMIall, file = "AnnV4_meta_PPMIall.csv", sep = ",", row.names=F,col.names=TRUE)
write.table(AnnV6_meta_PPMIall, file = "AnnV6_meta_PPMIall.csv", sep = ",", row.names=F,col.names=TRUE)
write.table(AnnV8_meta_PPMIall, file = "AnnV8_meta_PPMIall.csv", sep = ",", row.names=F,col.names=TRUE)





cmPD0.05 = metaPD0.05[metaPD0.05$MarkerName %in% commonDMP_PD_HC_0.05,]
cmHC0.05 = metaHC0.05[metaHC0.05$MarkerName %in% commonDMP_PD_HC_0.05,]
library(dplyr)
uniquePD0.05 <- anti_join(metaPD0.05, cmPD0.05)
uniqueHC0.05 <- anti_join(metaHC0.05, cmHC0.05)

write.table(uniquePD0.05, file = "uniquePD0.05.csv", sep = ",", row.names=FALSE,col.names=TRUE)
write.table(uniqueHC0.05, file = "uniqueHC0.05.csv", sep = ",", row.names=FALSE,col.names=TRUE)
write.table(cmPD0.05, file = "cmPD0.05.csv", sep = ",", row.names=FALSE,col.names=TRUE)
write.table(cmHC0.05, file = "cmHC0.05.csv", sep = ",", row.names=FALSE,col.names=TRUE)

#___________________________________________________________________________________
rm(metaHC)
metaPD2.4=read.table("2.4e-7PDMETAANALYSIS1.TBL",sep="\t",header=T) #4298
rownames(metaPD2.4)=metaPD2.4[,1]
metaHC2.4=read.table("2.4e-7HCMETAANALYSIS1.TBL",sep="\t",header=T) #1045
rownames(metaHC2.4)=metaHC2.4[,1]


metaPD2.4 <- read.csv("2.4metaPD.csv") #701
metaHC2.4 <- read.csv("2.4metaHC.csv") #297

commonDMP_PD_HC_2.4 <- intersect(metaPD2.4$MarkerName,metaHC2.4$MarkerName) #265

cmPD2.4 = metaPD2.4[metaPD2.4$MarkerName %in% commonDMP_PD_HC_2.4,]  #265
cmHC2.4 = metaHC2.4[metaHC2.4$MarkerName %in% commonDMP_PD_HC_2.4,]  #265
library(dplyr)
uniquePD2.4 <- anti_join(metaPD2.4, cmPD2.4)  #436
uniqueHC2.4 <- anti_join(metaHC2.4, cmHC2.4)  #32

write.table(uniquePD2.4, file = "uniquePD2.4.csv", sep = ",", row.names=FALSE,col.names=TRUE)
write.table(uniqueHC2.4, file = "uniqueHC2.4.csv", sep = ",", row.names=FALSE,col.names=TRUE)
write.table(cmPD2.4, file = "cmPD2.4.csv", sep = ",", row.names=FALSE,col.names=TRUE)
write.table(cmHC2.4, file = "cmHC2.4.csv", sep = ",", row.names=FALSE,col.names=TRUE)

#——————————————————————————————————————————————————————————————————————————————————————————————————————
metaPD0.05 <- read.csv("0.05metaPD.csv") #3025
metaHC0.05 <- read.csv("0.05metaHC.csv") #978
metaPD2.4 <- read.csv("2.4metaPD.csv") #701
metaHC2.4 <- read.csv("2.4metaHC.csv") #297

meta12PD <- intersect(metaPD0.05$MarkerName,metaPD2.4$MarkerName)
meta12HC <- intersect(metaHC0.05$MarkerName,metaHC2.4$MarkerName)

rm(metaHC)
metaPD0.05=read.table("0.05PDMETAANALYSIS1.TBL",sep="\t",header=T) #18644
rownames(metaPD0.05)=metaPD0.05[,1]
metaHC0.05=read.table("0.05HCMETAANALYSIS1.TBL",sep="\t",header=T) #4843
rownames(metaHC0.05)=metaHC0.05[,1]


metaPD0.05 <- read.csv("0.05metaPD.csv") #3024
metaHC0.05 <- read.csv("0.05metaHC.csv") #978

commonDMP_PD_HC_0.05 <- intersect(metaPD0.05$MarkerName,metaHC0.05$MarkerName) #825

cmPD0.05 = metaPD0.05[metaPD0.05$MarkerName %in% commonDMP_PD_HC_0.05,]  #825
cmHC0.05 = metaHC0.05[metaHC0.05$MarkerName %in% commonDMP_PD_HC_0.05,]  #825
library(dplyr)
uniquePD0.05 <- anti_join(metaPD0.05, cmPD0.05) #2199
uniqueHC0.05 <- anti_join(metaHC0.05, cmHC0.05) #153

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

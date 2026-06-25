### PDSex meta-analysis input preparation
PD331 <- read.csv("myDMPPD331_0.05_bacon.csv", row.names = 1)  #14045
PD813 <- read.csv("myDMPPD813_0.05_bacon.csv", row.names = 1)  #35631
HC236 <- read.csv("myDMPHC236_0.05_bacon.csv", row.names = 1)  #6173
HC770 <- read.csv("myDMPHC770_0.05_bacon.csv", row.names = 1)  #52468
PD331$CpGs <- row.names(PD331)
PD813$CpGs <- row.names(PD813)
HC236$CpGs <- row.names(HC236)
HC770$CpGs <- row.names(HC770)

PD331 <- subset(PD331, fdr.bacon < 0.05) #3504
PD813 <- subset(PD813, fdr.bacon < 0.05) #18165
HC236 <- subset(HC236, fdr.bacon < 0.05) #1244
HC770 <- subset(HC770, fdr.bacon < 0.05) #4593

write.table(PD331, file = "PD331_baconFDR0.05.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(PD813, file = "PD813_baconFDR0.05.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(HC236, file = "HC236_baconFDR0.05.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(HC770, file = "HC770_baconFDR0.05.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


PD331 <- read.csv("myDMPPD331_0.05_bacon.csv", row.names = 1)  #14045
PD813 <- read.csv("myDMPPD813_0.05_bacon.csv", row.names = 1)  #35631
HC236 <- read.csv("myDMPHC236_0.05_bacon.csv", row.names = 1)  #6173
HC770 <- read.csv("myDMPHC770_0.05_bacon.csv", row.names = 1)  #52468
PD331$CpGs <- row.names(PD331)
PD813$CpGs <- row.names(PD813)
HC236$CpGs <- row.names(HC236)
HC770$CpGs <- row.names(HC770)

PD331 <- subset(PD331, pValue.bacon < 2.4e-7) #745
PD813 <- subset(PD813, pValue.bacon < 2.4e-7) #4254
HC236 <- subset(HC236, pValue.bacon < 2.4e-7) #347
HC770 <- subset(HC770, pValue.bacon < 2.4e-7) #996

write.table(PD331, file = "PD331_baconpValue2.4e-7.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(PD813, file = "PD813_baconpValue2.4e-7.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(HC236, file = "HC236_baconpValue2.4e-7.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(HC770, file = "HC770_baconpValue2.4e-7.txt", sep = "\t", row.names = FALSE, col.names = TRUE)







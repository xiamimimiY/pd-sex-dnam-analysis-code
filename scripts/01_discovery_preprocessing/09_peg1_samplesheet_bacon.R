##
### ---------------
###
###
### ---------------


install.packages("openxlsx")
library(openxlsx)
a = read.xlsx("571samplesheet.xlsx", sheet = 1, startRow = 1, colNames = TRUE,
          skipEmptyRows = TRUE, rowNames = TRUE , detectDates = FALSE,
          rows = NULL, cols = NULL)



rm(phenoData)
## the probe design information of the array:RGChannelSet
manifest <- getManifest(rgSet)
manifest
head(getProbeInfo(manifest))



### ---------------------------------------------
library(bacon)
?bacon
y <- rnormmix(5000, c(0.9, 0.2, 1.3, 1, 4, 1))
class(y)   #[1] "numeric"
bc <- bacon(y)
bc
estimates(bc)  ##extract all estimated mixture parameters
inflation(bc)  #1.303211 ##extract only the inflation
bias(bc)   #0.1887784  ##extract only the bias
head(pval(bc))  ##extract bias and inflation corrected P-values
head(tstat(bc))  ##extract bias and inflation corrected test-statistics
traces(bc, burnin=FALSE)
posteriors(bc)
fit(bc, n=100)
plot(bc, type="hist")
plot(bc, type="qq")








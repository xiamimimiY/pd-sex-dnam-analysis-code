####
rm(list = ls())

library(MatrixEQTL)

base.dir = find.package('MatrixEQTL');


useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
snps_location_file_name = paste(base.dir, "/data/snpsloc.txt", sep="");

expression_file_name = paste(base.dir, "/data/GE.txt", sep="");
gene_location_file_name = paste(base.dir, "/data/geneloc.txt", sep="");

covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="");

output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 1e-2;

errorCovariance = numeric();
#errorCovariance = read.table("Sample_Data/errorCovariance.txt");

#Distance for local gene-SNP pairs
cisDist = 1e6;

snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 2000
snps$LoadFile( SNP_file_name )

gene = SlicedData$new();
gene$fileDelimiter = "\t";
gene$fileOmitCharacters = "NA"; 
gene$fileSkipRows = 1; 
gene$fileSkipColumns = 1; 
gene$fileSliceSize = 2000; 
gene$LoadFile(expression_file_name);

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t"; 
cvrt$fileOmitCharacters = "NA";
cvrt$fileSkipRows = 1; 
cvrt$fileSkipColumns = 1;
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name  = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = modelANOVA, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)
#Plot the Q-Q plot of local and distant p-values
plot(me)



###_____________________________________________________________________________________
###_____________________________________________________________________________________
library(MatrixEQTL)
rm(list = ls())
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 1e-2;

errorCovariance = numeric();

#Distance for local gene-cpg pairs
cisDist = 1e6;

DNAm_BL <- SlicedData$new()
DNAm_BL$fileDelimiter <- ","
DNAm_BL$fileOmitCharacters <- "NA"
DNAm_BL$fileSkipRows <- 1
DNAm_BL$fileSkipColumns = 1
DNAm_BL$fileSliceSize = 2000
DNAm_BL$LoadFile("DNAm_BL.csv")

eGene_BL <- SlicedData$new()
eGene_BL$fileDelimiter <- ","  
eGene_BL$fileOmitCharacters <- "NA" 
eGene_BL$fileSkipRows <- 1  
eGene_BL$fileSkipColumns = 1 
eGene_BL$fileSliceSize = 2000      
eGene_BL$LoadFile("eGene_BL.csv")

cova_BL = SlicedData$new();
cova_BL$fileDelimiter <- ","  
cova_BL$fileOmitCharacters = "NA";
cova_BL$fileSkipRows = 1; 
cova_BL$fileSkipColumns = 1;
cova_BL$LoadFile("cova_BL.csv")


cpgloc_BL <- read.csv("cpgloc_BL.csv", header = TRUE)
geneloc_BL <- read.csv("geneloc_BL.csv", header = TRUE)

me = Matrix_eQTL_main(
  snps = DNAm_BL, 
  gene = eGene_BL, 
  cvrt = cova_BL,
  output_file_name  = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = cpgloc_BL, 
  genepos = geneloc_BL,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);
class(me) #"list"       "MatrixEQTL"

unlink(output_file_name_tra);
unlink(output_file_name_cis);

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)
#Plot the Q-Q plot of local and distant p-values
plot(me)

write.csv(me$cis$eqtls, file = "cis_eQTMs_ageadjust.csv", row.names = FALSE)
write.csv(me$trans$eqtls, file = "trans_eQTMs_ageadjust.csv", row.names = FALSE)






##Creating a bed file as an input for the comb-p tool
##-------------------------------------------------------
colnames(uniquePDann)

dmr <- data.frame("chrom" = 	paste0("chr", uniquePDann$F_to_M.CHR),
                  "start" = 	uniquePDann[rownames(uniquePDann), "F_to_M.MAPINFO"],
                  "end" 	= 	uniquePDann[rownames(uniquePDann), "F_to_M.MAPINFO"]+1,
                  "pvalue" = 	uniquePDann$P.value)

dmr<-dmr[order(dmr[,1],dmr[,2]),]
write.table(dmr, file="forcombp.bed",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(dmr, file="combp.bed",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##-------------------------------------------------------
##Identification of Differentially methylated regions using the comb-p tool
##-------------------------------------------------------
# comb-p pipeline -c 4 --dist 500 --seed 1.0e-4 --anno hg19 -p  out meta.DMR.bed
# final command used for the manuscript:
# comb-p pipeline -c 4 --dist 750 --seed 0.05 --anno hg19 -p out1 combp.bed


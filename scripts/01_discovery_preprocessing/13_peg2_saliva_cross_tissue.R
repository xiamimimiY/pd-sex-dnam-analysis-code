DMP_PD <- read.csv("myDMPPD128_0.05_bacon.csv", header = TRUE, row.names = 1) #2674
DMP_HC <- read.csv("myDMPHC131_0.05_bacon.csv", header = TRUE, row.names = 1) #7414
DMP_HC <- DMP_HC[DMP_HC$fdr.bacon < 0.05, ] #362
colnames(DMP_PD)
common_rows <- intersect(rownames(DMP_PD), rownames(DMP_HC))
library(dplyr)
# Get the row names of DMP_PD and DMP_HC
row_names_PD <- rownames(DMP_PD)
row_names_HC <- rownames(DMP_HC)
# Find the rows unique to DMP_PD
DMP_PD_unique <- DMP_PD[!row_names_PD %in% row_names_HC, ]

# Set the file path
file_path <- "path/to/PDSex/meta & PPMI & saliva/uniquePD0.05.csv"
# Read the CSV file into a dataframe
uniquePDmeta <- read.csv(file_path, row.names = 1) #2199

common_rows1 <- intersect(rownames(DMP_PD_unique), rownames(uniquePDmeta))
class(common_rows1)
saliva491DMP <- DMP_PD_unique[rownames(DMP_PD_unique) %in% common_rows1, ]
# Get the column names of saliva491DMP
col_names <- colnames(saliva491DMP)
# Remove the "M_to_F." prefix from the column names
col_names <- gsub("M_to_F\\.", "", col_names)
# Set the modified column names back to saliva491DMP
colnames(saliva491DMP) <- col_names
# View the modified dataframe
print(saliva491DMP)

library(dplyr)
# Sort saliva491DMP by the fdr.bacon column in ascending order
saliva491DMP <- saliva491DMP[order(saliva491DMP$fdr.bacon), ]
DMP_PD_unique <- DMP_PD_unique[order(DMP_PD_unique$fdr.bacon), ]  
write.table(saliva491DMP, file = "saliva491uniqueDMP.csv", sep = ",") 
write.table(DMP_PD_unique, file = "saliva2373DMP_PD_unique.csv", sep = ",") 

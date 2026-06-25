
##2026-6-17
############################################################
## R2-9: Generate cell proportions and rerun MatrixEQTL eQTM
## Model: expression ~ DNAm + Age + estimated blood cell proportions
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

install.packages("MatrixEQTL")
library(MatrixEQTL)

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

library(dplyr)
library(openxlsx)

############################################################
## 1. Load ChAMP refbase result and inspect CellFraction
############################################################

load("ppmi140_champ_myRefbase_BL_PD213.Rdata")

cat("Objects loaded:\n")
print(ls())

cat("\nNames in myRefbasePD:\n")
print(names(myRefbasePD))

cat("\nDimension of CorrectedBeta:\n")
print(dim(myRefbasePD$CorrectedBeta))

cat("\nDimension of CellFraction:\n")
print(dim(myRefbasePD$CellFraction))

cat("\nHead of CellFraction:\n")
print(head(myRefbasePD$CellFraction))

cat("\nColumn names of pD:\n")
print(colnames(pD))

cat("\nHead of pD:\n")
print(head(pD[, 1:min(10, ncol(pD))]))

############################################################
## 2. Convert myRefbasePD$CellFraction to cellprop_BL.csv
############################################################

cell_raw <- myRefbasePD$CellFraction

cell_keywords <- c(
  "CD8T", "CD4T", "NK", "Bcell", "B_cell",
  "Mono", "Gran", "Neu", "Neutrophil", "Neutrophils",
  "Eos", "Bas"
)

if (sum(cell_keywords %in% colnames(cell_raw)) >= 3) {
  cell_df <- as.data.frame(cell_raw, check.names = FALSE)
  cell_df$CellID <- rownames(cell_df)
} else if (sum(cell_keywords %in% rownames(cell_raw)) >= 3) {
  cell_df <- as.data.frame(t(cell_raw), check.names = FALSE)
  cell_df$CellID <- rownames(cell_df)
} else {
  stop("Cannot identify cell-type columns/rows in myRefbasePD$CellFraction. Please inspect rownames/colnames.")
}

names(cell_df)[names(cell_df) == "B_cell"] <- "Bcell"
names(cell_df)[names(cell_df) == "Neutrophil"] <- "Neu"
names(cell_df)[names(cell_df) == "Neutrophils"] <- "Neu"
names(cell_df)[names(cell_df) == "Granulocyte"] <- "Gran"
names(cell_df)[names(cell_df) == "Granulocytes"] <- "Gran"

cell_cols <- intersect(
  c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "Neu", "Eos", "Bas"),
  names(cell_df)
)

cat("\nDetected cell fraction columns:\n")
print(cell_cols)

if (length(cell_cols) < 3) {
  stop("Too few cell proportion columns detected.")
}

############################################################
## 3. Match CellFraction samples to eQTM samples
############################################################

dnam_header <- read.csv("DNAm_BL.csv", nrows = 1, check.names = FALSE)
expr_header <- read.csv("eGene_BL.csv", nrows = 1, check.names = FALSE)

dnam_samples <- colnames(dnam_header)[-1]
expr_samples <- colnames(expr_header)[-1]

cat("\nNumber of DNAm samples:", length(dnam_samples), "\n")
cat("Number of expression samples:", length(expr_samples), "\n")

if (!identical(dnam_samples, expr_samples)) {
  stop("DNAm_BL.csv and eGene_BL.csv sample orders are not identical.")
}

eqtm_samples <- dnam_samples

direct_match_n <- sum(eqtm_samples %in% cell_df$CellID)
cat("\nDirect match between eQTM samples and CellFraction CellID:", direct_match_n, "/", length(eqtm_samples), "\n")

if (direct_match_n == length(eqtm_samples)) {
  
  cellprop_BL <- cell_df[match(eqtm_samples, cell_df$CellID), c("CellID", cell_cols), drop = FALSE]
  colnames(cellprop_BL)[1] <- "SampleID"
  
} else {
  
  pD_df <- as.data.frame(pD, check.names = FALSE)
  pD_df$CellID <- rownames(pD_df)
  
  pD_df <- pD_df[pD_df$CellID %in% cell_df$CellID, , drop = FALSE]
  
  candidate_id_cols <- colnames(pD_df)
  
  match_summary <- data.frame(
    column = candidate_id_cols,
    matched_n = sapply(candidate_id_cols, function(cc) {
      sum(eqtm_samples %in% as.character(pD_df[[cc]]))
    }),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(matched_n))
  
  cat("\nPotential sample-ID columns in pD:\n")
  print(head(match_summary, 20))
  
  best_col <- match_summary$column[1]
  best_n <- match_summary$matched_n[1]
  
  if (best_n != length(eqtm_samples)) {
    cat("\nFirst eQTM sample IDs:\n")
    print(head(eqtm_samples, 20))
    cat("\nFirst CellFraction CellIDs:\n")
    print(head(cell_df$CellID, 20))
    cat("\nTop pD matching columns:\n")
    print(head(match_summary, 20))
    stop("Cannot fully match eQTM samples to CellFraction. Please check sample IDs.")
  }
  
  cat("\nUsing pD column for sample matching:", best_col, "\n")
  
  map_df <- pD_df[, c("CellID", best_col), drop = FALSE]
  colnames(map_df) <- c("CellID", "SampleID")
  map_df$SampleID <- as.character(map_df$SampleID)
  
  cell_mapped <- left_join(map_df, cell_df, by = "CellID")
  
  cellprop_BL <- cell_mapped[match(eqtm_samples, cell_mapped$SampleID), c("SampleID", cell_cols), drop = FALSE]
}

cat("\nFinal cellprop_BL dimension:\n")
print(dim(cellprop_BL))

cat("\nHead of cellprop_BL:\n")
print(head(cellprop_BL))

if (!identical(cellprop_BL$SampleID, eqtm_samples)) {
  stop("cellprop_BL sample order does not match DNAm_BL.csv/eGene_BL.csv.")
}

if (any(is.na(cellprop_BL))) {
  stop("NA detected in cellprop_BL. Please check sample matching or cell fraction values.")
}

write.csv(cellprop_BL, "cellprop_BL.csv", row.names = FALSE)
cat("\nSaved: cellprop_BL.csv\n")

############################################################
## 4. Build new MatrixEQTL covariate file: Age + cell fractions
############################################################

cova_old <- read.csv("cova_BL.csv", check.names = FALSE)

cat("\nOriginal cova_BL.csv:\n")
print(cova_old[, 1:min(6, ncol(cova_old))])

cov_name_col <- colnames(cova_old)[1]
cov_names <- as.character(cova_old[[1]])

age_row <- which(tolower(cov_names) == "age")

if (length(age_row) != 1) {
  stop("Cannot uniquely identify Age row in cova_BL.csv.")
}

age_values <- as.numeric(cova_old[age_row, -1])
names(age_values) <- colnames(cova_old)[-1]

if (!all(eqtm_samples %in% names(age_values))) {
  stop("Some eQTM samples are missing from cova_BL.csv.")
}

age_values <- age_values[eqtm_samples]

cell_cols_model <- cell_cols

if ("Gran" %in% cell_cols_model) {
  cell_cols_model <- setdiff(cell_cols_model, "Gran")
} else if ("Neu" %in% cell_cols_model) {
  cell_cols_model <- setdiff(cell_cols_model, "Neu")
} else {
  cell_cols_model <- cell_cols_model[-length(cell_cols_model)]
}

cat("\nCell fractions used in MatrixEQTL covariates:\n")
print(cell_cols_model)

cell_mat <- as.matrix(cellprop_BL[, cell_cols_model, drop = FALSE])
mode(cell_mat) <- "numeric"
rownames(cell_mat) <- cellprop_BL$SampleID

cov_mat <- rbind(
  Age = age_values,
  t(cell_mat)
)

colnames(cov_mat) <- eqtm_samples

if (any(is.na(cov_mat))) {
  stop("NA detected in final covariate matrix.")
}

cova_age_cell <- data.frame(
  Covariate = rownames(cov_mat),
  cov_mat,
  check.names = FALSE
)

write.csv(cova_age_cell, "cova_BL_age_cell.csv", row.names = FALSE, quote = FALSE)

cat("\nSaved: cova_BL_age_cell.csv\n")
cat("Final covariates:\n")
print(rownames(cov_mat))

############################################################
## 5. Run MatrixEQTL: Age + cell fractions
############################################################

useModel <- modelLINEAR
errorCovariance <- numeric()
cisDist <- 1e6

pvOutputThreshold_cis <- 1

pvOutputThreshold_tra <- 0

output_file_name_cis <- "cis_eQTMs_age_cell_adjusted_raw.txt"
output_file_name_tra <- tempfile()

DNAm_BL <- SlicedData$new()
DNAm_BL$fileDelimiter <- ","
DNAm_BL$fileOmitCharacters <- "NA"
DNAm_BL$fileSkipRows <- 1
DNAm_BL$fileSkipColumns <- 1
DNAm_BL$fileSliceSize <- 2000
DNAm_BL$LoadFile("DNAm_BL.csv")

eGene_BL <- SlicedData$new()
eGene_BL$fileDelimiter <- ","
eGene_BL$fileOmitCharacters <- "NA"
eGene_BL$fileSkipRows <- 1
eGene_BL$fileSkipColumns <- 1
eGene_BL$fileSliceSize <- 2000
eGene_BL$LoadFile("eGene_BL.csv")

cova_BL_age_cell <- SlicedData$new()
cova_BL_age_cell$fileDelimiter <- ","
cova_BL_age_cell$fileOmitCharacters <- "NA"
cova_BL_age_cell$fileSkipRows <- 1
cova_BL_age_cell$fileSkipColumns <- 1
cova_BL_age_cell$fileSliceSize <- 2000
cova_BL_age_cell$LoadFile("cova_BL_age_cell.csv")

cpgloc_BL <- read.csv("cpgloc_BL.csv", header = TRUE, stringsAsFactors = FALSE)
geneloc_BL <- read.csv("geneloc_BL.csv", header = TRUE, stringsAsFactors = FALSE)

cat("\nCpG location columns:\n")
print(names(cpgloc_BL))
cat("\nGene location columns:\n")
print(names(geneloc_BL))

if (!all(c("snpid", "chr", "pos") %in% names(cpgloc_BL))) {
  stop("cpgloc_BL.csv must contain columns: snpid, chr, pos")
}

if (!all(c("geneid", "chr", "s1", "s2") %in% names(geneloc_BL))) {
  stop("geneloc_BL.csv must contain columns: geneid, chr, s1, s2")
}

me_age_cell <- Matrix_eQTL_main(
  snps = DNAm_BL,
  gene = eGene_BL,
  cvrt = cova_BL_age_cell,
  output_file_name = output_file_name_tra,
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
  noFDRsaveMemory = FALSE
)

unlink(output_file_name_tra)

cat("\nMatrixEQTL completed in:", me_age_cell$time.in.sec, "seconds\n")

############################################################
## 6. Save cis-eQTM results
############################################################

cis_age_cell <- me_age_cell$cis$eqtls

cat("\nColumns in cis_age_cell:\n")
print(names(cis_age_cell))

if (!"FDR" %in% names(cis_age_cell)) {
  cis_age_cell$FDR <- p.adjust(cis_age_cell$pvalue, method = "BH")
}

cis_age_cell <- cis_age_cell[order(cis_age_cell$FDR, cis_age_cell$pvalue), ]

write.csv(cis_age_cell, "cis_eQTMs_age_cell_adjusted_all.csv", row.names = FALSE)

cis_age_cell_sig <- cis_age_cell %>%
  filter(FDR < 0.05)

write.csv(cis_age_cell_sig, "cis_eQTMs_age_cell_adjusted_FDR005.csv", row.names = FALSE)

cat("\nAge + cell adjusted cis-eQTMs FDR<0.05:", nrow(cis_age_cell_sig), "\n")

############################################################
## 7. Extract and compare ZNF727 six cis-eQTMs
############################################################

znf727_cpgs <- c(
  "cg06098368",
  "cg03124146",
  "cg12479444",
  "cg20067334",
  "cg14285533",
  "cg26911611"
)

znf727_gene_candidates <- c("ZNF727", "ENSG00000214652")

znf727_age_cell <- cis_age_cell %>%
  filter(snps %in% znf727_cpgs, gene %in% znf727_gene_candidates) %>%
  arrange(snps)

write.csv(znf727_age_cell, "ZNF727_six_cis_eQTMs_age_cell_adjusted.csv", row.names = FALSE)

cat("\nZNF727 cis-eQTMs found after age + cell adjustment:", nrow(znf727_age_cell), "\n")
cat("ZNF727 cis-eQTMs with FDR<0.05 after age + cell adjustment:",
    sum(znf727_age_cell$FDR < 0.05, na.rm = TRUE), "\n")

print(znf727_age_cell)

############################################################
## 8. Compare with original age-adjusted eQTM result
############################################################

if (file.exists("cis_eQTMs_ageadjust.csv")) {
  
  cis_age_only <- read.csv("cis_eQTMs_ageadjust.csv", check.names = FALSE)
  
  if (!"FDR" %in% names(cis_age_only)) {
    cis_age_only$FDR <- p.adjust(cis_age_only$pvalue, method = "BH")
  }
  
  znf727_age_only <- cis_age_only %>%
    filter(snps %in% znf727_cpgs, gene %in% znf727_gene_candidates)
  
  common_cols <- intersect(c("snps", "gene", "statistic", "pvalue", "FDR", "beta"), names(znf727_age_only))
  
  znf727_age_only2 <- znf727_age_only[, common_cols, drop = FALSE]
  znf727_age_cell2 <- znf727_age_cell[, common_cols, drop = FALSE]
  
  names(znf727_age_only2) <- paste0(names(znf727_age_only2), "_age_only")
  names(znf727_age_cell2) <- paste0(names(znf727_age_cell2), "_age_cell")
  
  names(znf727_age_only2)[names(znf727_age_only2) == "snps_age_only"] <- "snps"
  names(znf727_age_only2)[names(znf727_age_only2) == "gene_age_only"] <- "gene"
  names(znf727_age_cell2)[names(znf727_age_cell2) == "snps_age_cell"] <- "snps"
  names(znf727_age_cell2)[names(znf727_age_cell2) == "gene_age_cell"] <- "gene"
  
  znf727_compare <- full_join(
    znf727_age_only2,
    znf727_age_cell2,
    by = c("snps", "gene")
  ) %>%
    mutate(
      Direction_consistent = sign(beta_age_only) == sign(beta_age_cell),
      Significant_age_only = FDR_age_only < 0.05,
      Significant_age_cell = FDR_age_cell < 0.05
    )
  
  write.csv(
    znf727_compare,
    "ZNF727_six_cis_eQTMs_age_only_vs_age_cell.csv",
    row.names = FALSE
  )
  
  cat("\nZNF727 age-only vs age+cell comparison:\n")
  print(znf727_compare)
}

############################################################
## 9. Summary for R2-9 response
############################################################

summary_counts <- data.frame(
  Analysis = c(
    "cis_eQTMs_FDR005_age_cell",
    "ZNF727_six_cis_eQTMs_expected",
    "ZNF727_six_cis_eQTMs_found_age_cell",
    "ZNF727_six_cis_eQTMs_FDR005_age_cell"
  ),
  N = c(
    nrow(cis_age_cell_sig),
    length(znf727_cpgs),
    nrow(znf727_age_cell),
    sum(znf727_age_cell$FDR < 0.05, na.rm = TRUE)
  )
)

print(summary_counts)

write.csv(summary_counts, "R2_9_eQTM_age_cell_summary_counts.csv", row.names = FALSE)

openxlsx::write.xlsx(
  list(
    Summary_counts = summary_counts,
    ZNF727_age_cell = znf727_age_cell,
    ZNF727_age_only_vs_age_cell = if (exists("znf727_compare")) znf727_compare else data.frame(),
    cis_age_cell_FDR005 = cis_age_cell_sig,
    cis_age_cell_all = cis_age_cell
  ),
  file = "R2_9_eQTM_age_cell_adjusted_summary.xlsx",
  overwrite = TRUE
)

cat("\nR2-9 analysis completed.\n")
cat("Key files generated:\n")
cat("1. cellprop_BL.csv\n")
cat("2. cova_BL_age_cell.csv\n")
cat("3. cis_eQTMs_age_cell_adjusted_all.csv\n")
cat("4. cis_eQTMs_age_cell_adjusted_FDR005.csv\n")
cat("5. ZNF727_six_cis_eQTMs_age_cell_adjusted.csv\n")
cat("6. ZNF727_six_cis_eQTMs_age_only_vs_age_cell.csv\n")
cat("7. R2_9_eQTM_age_cell_summary_counts.csv\n")
cat("8. R2_9_eQTM_age_cell_adjusted_summary.xlsx\n")



############################################################
## R2-9: Generate ZNF727 comparison table only
## No need to rerun MatrixEQTL
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

############################################################
## 1. Read age-only and age+cell-adjusted cis-eQTM results
############################################################

cis_age_only <- read.csv("cis_eQTMs_ageadjust.csv", check.names = FALSE)
cis_age_cell <- read.csv("cis_eQTMs_age_cell_adjusted_all.csv", check.names = FALSE)

## Remove empty or NA column names
clean_empty_names <- function(df) {
  nm <- names(df)
  keep <- !(is.na(nm) | nm == "")
  df <- df[, keep, drop = FALSE]
  names(df) <- gsub("\\.", "_", names(df))
  return(df)
}

cis_age_only <- clean_empty_names(cis_age_only)
cis_age_cell <- clean_empty_names(cis_age_cell)

## Standardize possible column names
standardize_eqtm_names <- function(df) {
  names(df) <- gsub("t_stat", "statistic", names(df))
  names(df) <- gsub("P_value", "pvalue", names(df))
  names(df) <- gsub("p_value", "pvalue", names(df))
  names(df) <- gsub("fdr", "FDR", names(df), ignore.case = FALSE)
  return(df)
}

cis_age_only <- standardize_eqtm_names(cis_age_only)
cis_age_cell <- standardize_eqtm_names(cis_age_cell)

## Check required columns
required_cols <- c("snps", "gene", "statistic", "pvalue", "FDR", "beta")

cat("Columns in cis_age_only:\n")
print(names(cis_age_only))

cat("Columns in cis_age_cell:\n")
print(names(cis_age_cell))

if (!all(required_cols %in% names(cis_age_only))) {
  stop("cis_eQTMs_ageadjust.csv is missing required columns: ",
       paste(setdiff(required_cols, names(cis_age_only)), collapse = ", "))
}

if (!all(required_cols %in% names(cis_age_cell))) {
  stop("cis_eQTMs_age_cell_adjusted_all.csv is missing required columns: ",
       paste(setdiff(required_cols, names(cis_age_cell)), collapse = ", "))
}

############################################################
## 2. Extract six ZNF727 cis-eQTMs
############################################################

znf727_cpgs <- c(
  "cg06098368",
  "cg03124146",
  "cg12479444",
  "cg20067334",
  "cg14285533",
  "cg26911611"
)

znf727_gene_candidates <- c("ZNF727", "ENSG00000214652")

age_only <- cis_age_only[
  cis_age_only$snps %in% znf727_cpgs &
    cis_age_only$gene %in% znf727_gene_candidates,
  required_cols
]

age_cell <- cis_age_cell[
  cis_age_cell$snps %in% znf727_cpgs &
    cis_age_cell$gene %in% znf727_gene_candidates,
  required_cols
]

## Rename columns for comparison
names(age_only)[3:6] <- paste0(names(age_only)[3:6], "_age_only")
names(age_cell)[3:6] <- paste0(names(age_cell)[3:6], "_age_cell")

znf727_compare <- merge(
  age_only,
  age_cell,
  by = c("snps", "gene"),
  all = TRUE
)

## Keep the expected CpG order
znf727_compare <- znf727_compare[
  match(znf727_cpgs, znf727_compare$snps),
]

znf727_compare$Direction_consistent <- sign(znf727_compare$beta_age_only) == sign(znf727_compare$beta_age_cell)
znf727_compare$Significant_age_only <- znf727_compare$FDR_age_only < 0.05
znf727_compare$Significant_age_cell <- znf727_compare$FDR_age_cell < 0.05

############################################################
## 3. Summary table
############################################################

cis_age_cell_sig <- cis_age_cell[cis_age_cell$FDR < 0.05, ]

summary_counts <- data.frame(
  Analysis = c(
    "cis_eQTMs_FDR005_age_cell",
    "ZNF727_six_cis_eQTMs_expected",
    "ZNF727_six_cis_eQTMs_found_age_cell",
    "ZNF727_six_cis_eQTMs_FDR005_age_cell",
    "ZNF727_six_cis_eQTMs_direction_consistent"
  ),
  N = c(
    nrow(cis_age_cell_sig),
    length(znf727_cpgs),
    sum(!is.na(znf727_compare$FDR_age_cell)),
    sum(znf727_compare$FDR_age_cell < 0.05, na.rm = TRUE),
    sum(znf727_compare$Direction_consistent, na.rm = TRUE)
  )
)

############################################################
## 4. Save outputs
############################################################

write.csv(
  znf727_compare,
  "ZNF727_six_cis_eQTMs_age_only_vs_age_cell.csv",
  row.names = FALSE
)

write.csv(
  summary_counts,
  "R2_9_eQTM_age_cell_summary_counts_final.csv",
  row.names = FALSE
)

openxlsx::write.xlsx(
  list(
    Summary_counts = summary_counts,
    ZNF727_age_only_vs_age_cell = znf727_compare,
    cis_age_cell_FDR005 = cis_age_cell_sig
  ),
  file = "Supplementary_eTable14_eQTM_age_cell_adjusted_sensitivity.xlsx",
  overwrite = TRUE
)

############################################################
## 5. Print for copy
############################################################

cat("\nSummary counts:\n")
print(summary_counts)

cat("\nZNF727 age-only vs age+cell comparison:\n")
print(znf727_compare)

cat("\nFiles generated:\n")
cat("1. ZNF727_six_cis_eQTMs_age_only_vs_age_cell.csv\n")
cat("2. R2_9_eQTM_age_cell_summary_counts_final.csv\n")
cat("3. Supplementary_eTable14_eQTM_age_cell_adjusted_sensitivity.xlsx\n")


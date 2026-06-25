
###————————————2026-6-14————————————————###
############################################################
## R2-1: Background comparison of DMP sets
## Purpose:
## Compare PD-unique DMPs, PD-HC common DMPs, and HC-unique DMPs
## Analyses:
## 1. Set counts
## 2. Direction distribution
## 3. Gene feature distribution
## 4. CpG feature distribution
## 5. Chi-square tests for feature distributions
## 6. GO enrichment using missMethyl::gometh
## 7. Export summary Excel
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## R2-1: Background comparison of DMP sets
## Updated version: no count() used
## Input files:
## 1. 2199_PD-unique-DMPs.xlsx
## 2. 978_HC-sex-DMPs.xlsx
## 3. 826_PD-HC-cmDMPs.xlsx
############################################################
############################################################
## 0. Working directory
############################################################
# setwd("path/to/files")

############################################################
## 1. Load / install packages
############################################################

cran_pkgs <- c("readxl", "dplyr", "stringr", "openxlsx", "ggplot2")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(readxl)
library(dplyr)
library(stringr)
library(openxlsx)
library(ggplot2)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "missMethyl",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "minfi"
)

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)

############################################################
## 2. Read cleaned Excel files
############################################################
## Required Excel columns:
## CpG, Effect, StdErr, P.value, Direction, CHR, Position,
## Gene, Gene feature, CpG features

read_dmp_clean <- function(file, set_name) {
  
  df <- readxl::read_excel(file)
  
  message("Reading file: ", file)
  message("Original column names:")
  print(names(df))
  
  names(df) <- names(df) |>
    stringr::str_replace_all("\\.", "_") |>
    stringr::str_replace_all(" ", "_")
  
  message("Standardized column names:")
  print(names(df))
  
  required_cols <- c(
    "CpG", "Effect", "StdErr", "P_value", "Direction",
    "CHR", "Position", "Gene", "Gene_feature", "CpG_features"
  )
  
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", file, ": ",
      paste(missing_cols, collapse = ", "),
      "\nPlease check Excel column names."
    )
  }
  
  df_clean <- df |>
    dplyr::transmute(
      CpG = as.character(CpG),
      Effect = as.numeric(Effect),
      StdErr = as.numeric(StdErr),
      P_value = as.numeric(P_value),
      Direction = as.character(Direction),
      CHR = as.character(CHR),
      Position = as.numeric(Position),
      Gene = as.character(Gene),
      Gene_feature = as.character(Gene_feature),
      CpG_feature = as.character(CpG_features),
      Set = set_name
    ) |>
    dplyr::filter(!is.na(CpG), CpG != "") |>
    dplyr::mutate(
      Gene = ifelse(is.na(Gene) | Gene == "" | Gene == "NA", NA, Gene),
      Gene_feature = ifelse(
        is.na(Gene_feature) | Gene_feature == "" | Gene_feature == "NA",
        "IGR",
        Gene_feature
      ),
      CpG_feature = ifelse(
        is.na(CpG_feature) | CpG_feature == "" | CpG_feature == "NA",
        "unknown",
        CpG_feature
      )
    ) |>
    dplyr::distinct(CpG, .keep_all = TRUE)
  
  message(file, " read successfully: ", nrow(df_clean), " unique CpGs.")
  return(df_clean)
}

pd_unique <- read_dmp_clean("2199_PD-unique-DMPs.xlsx", "PD_unique")
hc_all <- read_dmp_clean("978_HC-sex-DMPs.xlsx", "HC_all")
pd_hc_common <- read_dmp_clean("826_PD-HC-cmDMPs.xlsx", "PD_HC_common")

############################################################
## 3. Generate true HC-unique DMPs
############################################################
## HC-unique = all HC sex-DMPs minus PD-HC common DMPs

hc_unique <- hc_all |>
  dplyr::filter(!CpG %in% pd_hc_common$CpG) |>
  dplyr::mutate(Set = "HC_unique")

set_counts <- data.frame(
  Set = c("PD_unique", "PD_HC_common", "HC_all", "HC_unique"),
  N = c(
    nrow(pd_unique),
    nrow(pd_hc_common),
    nrow(hc_all),
    nrow(hc_unique)
  )
)

print(set_counts)

write.csv(set_counts, "R2_1_DMP_set_counts.csv", row.names = FALSE)

message("Expected counts:")
message("PD_unique = 2199")
message("PD_HC_common = 826")
message("HC_all = 978")
message("HC_unique = 152")

############################################################
## 4. Combine three sets for main comparison
############################################################

dmp3 <- dplyr::bind_rows(
  pd_unique |> dplyr::mutate(Set = "PD_unique"),
  pd_hc_common |> dplyr::mutate(Set = "PD_HC_common"),
  hc_unique |> dplyr::mutate(Set = "HC_unique")
)

write.csv(dmp3, "R2_1_DMP_three_sets_master.csv", row.names = FALSE)

############################################################
## 5. Direction distribution
############################################################

direction_summary <- dmp3 |>
  dplyr::mutate(
    Direction_class = dplyr::case_when(
      Effect > 0 ~ "Male_hypermethylated",
      Effect < 0 ~ "Male_hypomethylated",
      TRUE ~ "No_direction"
    )
  ) |>
  dplyr::group_by(Set, Direction_class) |>
  dplyr::summarise(N = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(Set) |>
  dplyr::mutate(Percent = N / sum(N) * 100) |>
  dplyr::ungroup() |>
  dplyr::arrange(Set, Direction_class)

print(direction_summary)

write.csv(direction_summary, "R2_1_DMP_direction_summary.csv", row.names = FALSE)

############################################################
## 6. Simplify gene feature and CpG feature annotations
############################################################

simplify_gene_feature <- function(x) {
  
  if (is.na(x) || x == "" || toupper(x) == "NA") {
    return("Intergenic")
  }
  
  x0 <- tolower(gsub("'", "", x))
  vals <- trimws(unlist(strsplit(x0, ";|,|\\|")))
  
  if (any(vals %in% c("tss200", "tss1500", "5utr", "1stexon"))) {
    return("Promoter")
  }
  
  if (any(vals %in% c("body"))) {
    return("Gene_body")
  }
  
  if (any(vals %in% c("3utr"))) {
    return("3UTR")
  }
  
  if (any(vals %in% c("igr", "intergenic"))) {
    return("Intergenic")
  }
  
  return("Other")
}

simplify_cpg_feature <- function(x) {
  
  if (is.na(x) || x == "" || toupper(x) == "NA") {
    return("Unknown")
  }
  
  x0 <- tolower(x)
  
  if (stringr::str_detect(x0, "island")) {
    return("Island")
  }
  
  if (stringr::str_detect(x0, "shore")) {
    return("Shore")
  }
  
  if (stringr::str_detect(x0, "shelf")) {
    return("Shelf")
  }
  
  if (stringr::str_detect(x0, "open")) {
    return("OpenSea")
  }
  
  return("Unknown")
}

dmp3 <- dmp3 |>
  dplyr::mutate(
    Gene_feature_simple = sapply(Gene_feature, simplify_gene_feature),
    CpG_feature_simple = sapply(CpG_feature, simplify_cpg_feature)
  )

write.csv(dmp3, "R2_1_DMP_three_sets_master_with_clean_features.csv", row.names = FALSE)

############################################################
## 7. Gene feature and CpG feature distribution
############################################################

gene_feature_dist <- dmp3 |>
  dplyr::group_by(Set, Gene_feature_simple) |>
  dplyr::summarise(N = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(Set) |>
  dplyr::mutate(Percent = N / sum(N) * 100) |>
  dplyr::ungroup() |>
  dplyr::arrange(Set, Gene_feature_simple)

cpg_feature_dist <- dmp3 |>
  dplyr::group_by(Set, CpG_feature_simple) |>
  dplyr::summarise(N = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(Set) |>
  dplyr::mutate(Percent = N / sum(N) * 100) |>
  dplyr::ungroup() |>
  dplyr::arrange(Set, CpG_feature_simple)

print(gene_feature_dist)
print(cpg_feature_dist)

write.csv(gene_feature_dist, "R2_1_gene_feature_distribution.csv", row.names = FALSE)
write.csv(cpg_feature_dist, "R2_1_CpG_feature_distribution.csv", row.names = FALSE)

############################################################
## 8. Chi-square tests for feature distribution
############################################################

gene_feature_mat <- table(dmp3$Set, dmp3$Gene_feature_simple)
cpg_feature_mat <- table(dmp3$Set, dmp3$CpG_feature_simple)

gene_feature_chisq <- chisq.test(gene_feature_mat)
cpg_feature_chisq <- chisq.test(cpg_feature_mat)

print(gene_feature_mat)
print(gene_feature_chisq)

print(cpg_feature_mat)
print(cpg_feature_chisq)

capture.output(
  list(
    gene_feature_table = gene_feature_mat,
    gene_feature_chisq = gene_feature_chisq,
    cpg_feature_table = cpg_feature_mat,
    cpg_feature_chisq = cpg_feature_chisq
  ),
  file = "R2_1_feature_chisq_tests.txt"
)

############################################################
## 9. Per-feature Fisher enrichment tests
############################################################
## Each set vs the other two sets combined.

feature_fisher <- function(df, feature_col) {
  
  sets <- unique(df$Set)
  features <- unique(df[[feature_col]])
  
  res <- list()
  k <- 1
  
  for (s in sets) {
    for (f in features) {
      
      in_set <- df$Set == s
      has_feature <- df[[feature_col]] == f
      
      a <- sum(in_set & has_feature)
      b <- sum(in_set & !has_feature)
      c <- sum(!in_set & has_feature)
      d <- sum(!in_set & !has_feature)
      
      ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
      
      res[[k]] <- data.frame(
        Set = s,
        Feature = f,
        In_set_with_feature = a,
        In_set_without_feature = b,
        Other_sets_with_feature = c,
        Other_sets_without_feature = d,
        OR = unname(ft$estimate),
        P_value = ft$p.value
      )
      
      k <- k + 1
    }
  }
  
  out <- dplyr::bind_rows(res)
  out$FDR <- p.adjust(out$P_value, method = "BH")
  out <- out |> dplyr::arrange(FDR, P_value)
  
  return(out)
}

gene_feature_fisher <- feature_fisher(dmp3, "Gene_feature_simple")
cpg_feature_fisher <- feature_fisher(dmp3, "CpG_feature_simple")

write.csv(gene_feature_fisher, "R2_1_gene_feature_fisher_enrichment.csv", row.names = FALSE)
write.csv(cpg_feature_fisher, "R2_1_CpG_feature_fisher_enrichment.csv", row.names = FALSE)

############################################################
## 10. GO enrichment using missMethyl::gometh
############################################################
## gometh accounts for CpG-number bias across genes.

ann <- as.data.frame(minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
ann$CpG <- rownames(ann)

## Background option:
## If you have the actual tested CpGs, create:
## background_tested_CpGs.csv
## with one column named CpG.
## Otherwise, all autosomal 450K probes are used.

if (file.exists("background_tested_CpGs.csv")) {
  
  background_df <- read.csv("background_tested_CpGs.csv")
  background <- unique(as.character(background_df$CpG))
  background <- intersect(background, ann$CpG)
  message("Using actual tested CpGs as background: ", length(background))
  
} else {
  
  background <- ann$CpG[!ann$chr %in% c("chrX", "chrY")]
  message("Using 450K autosomal CpGs as background: ", length(background))
  
}

run_gometh <- function(cpgs, set_name) {
  
  sig <- unique(as.character(cpgs))
  sig <- intersect(sig, background)
  
  message(set_name, ": ", length(sig), " CpGs used for GO enrichment.")
  
  res <- missMethyl::gometh(
    sig.cpg = sig,
    all.cpg = background,
    collection = "GO",
    array.type = "450K",
    plot.bias = FALSE
  )
  
  res$Set <- set_name
  res <- res |> dplyr::arrange(FDR, P.DE)
  
  return(res)
}

go_pd_unique <- run_gometh(pd_unique$CpG, "PD_unique")
go_pd_hc_common <- run_gometh(pd_hc_common$CpG, "PD_HC_common")
go_hc_unique <- run_gometh(hc_unique$CpG, "HC_unique")

write.csv(go_pd_unique, "R2_1_GO_missMethyl_PD_unique_DMPs.csv", row.names = FALSE)
write.csv(go_pd_hc_common, "R2_1_GO_missMethyl_PD_HC_common_DMPs.csv", row.names = FALSE)
write.csv(go_hc_unique, "R2_1_GO_missMethyl_HC_unique_DMPs.csv", row.names = FALSE)

go_counts <- data.frame(
  Set = c("PD_unique", "PD_HC_common", "HC_unique"),
  N_CpGs = c(
    length(intersect(pd_unique$CpG, background)),
    length(intersect(pd_hc_common$CpG, background)),
    length(intersect(hc_unique$CpG, background))
  ),
  GO_terms_FDR_0_05 = c(
    sum(go_pd_unique$FDR < 0.05, na.rm = TRUE),
    sum(go_pd_hc_common$FDR < 0.05, na.rm = TRUE),
    sum(go_hc_unique$FDR < 0.05, na.rm = TRUE)
  ),
  GO_terms_FDR_0_10 = c(
    sum(go_pd_unique$FDR < 0.10, na.rm = TRUE),
    sum(go_pd_hc_common$FDR < 0.10, na.rm = TRUE),
    sum(go_hc_unique$FDR < 0.10, na.rm = TRUE)
  )
)

print(go_counts)

write.csv(go_counts, "R2_1_GO_enrichment_counts.csv", row.names = FALSE)

go_top20_fdr005 <- dplyr::bind_rows(
  head(go_pd_unique |> dplyr::filter(FDR < 0.05), 20),
  head(go_pd_hc_common |> dplyr::filter(FDR < 0.05), 20),
  head(go_hc_unique |> dplyr::filter(FDR < 0.05), 20)
)

write.csv(go_top20_fdr005, "R2_1_GO_top20_FDR005_all_sets.csv", row.names = FALSE)

go_top20_nominal <- dplyr::bind_rows(
  head(go_pd_unique, 20),
  head(go_pd_hc_common, 20),
  head(go_hc_unique, 20)
)

write.csv(go_top20_nominal, "R2_1_GO_top20_nominal_all_sets.csv", row.names = FALSE)

############################################################
## 11. Optional plots
############################################################

p_gene <- ggplot2::ggplot(
  gene_feature_dist,
  ggplot2::aes(x = Set, y = Percent, fill = Gene_feature_simple)
) +
  ggplot2::geom_col(position = "stack") +
  ggplot2::theme_bw() +
  ggplot2::labs(
    x = NULL,
    y = "Percentage of DMPs (%)",
    fill = "Gene feature",
    title = "Gene feature distribution across DMP sets"
  ) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

ggplot2::ggsave("R2_1_gene_feature_distribution.pdf", p_gene, width = 7, height = 5)
ggplot2::ggsave("R2_1_gene_feature_distribution.png", p_gene, width = 7, height = 5, dpi = 300)

p_cpg <- ggplot2::ggplot(
  cpg_feature_dist,
  ggplot2::aes(x = Set, y = Percent, fill = CpG_feature_simple)
) +
  ggplot2::geom_col(position = "stack") +
  ggplot2::theme_bw() +
  ggplot2::labs(
    x = NULL,
    y = "Percentage of DMPs (%)",
    fill = "CpG feature",
    title = "CpG feature distribution across DMP sets"
  ) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

ggplot2::ggsave("R2_1_CpG_feature_distribution.pdf", p_cpg, width = 7, height = 5)
ggplot2::ggsave("R2_1_CpG_feature_distribution.png", p_cpg, width = 7, height = 5, dpi = 300)

############################################################
## 12. Save summary workbook
############################################################

openxlsx::write.xlsx(
  list(
    Set_counts = set_counts,
    Direction_summary = direction_summary,
    Gene_feature_distribution = gene_feature_dist,
    CpG_feature_distribution = cpg_feature_dist,
    Gene_feature_chisq_table = as.data.frame.matrix(gene_feature_mat),
    CpG_feature_chisq_table = as.data.frame.matrix(cpg_feature_mat),
    Gene_feature_Fisher = gene_feature_fisher,
    CpG_feature_Fisher = cpg_feature_fisher,
    GO_counts = go_counts,
    GO_top20_FDR005 = go_top20_fdr005,
    GO_top20_nominal = go_top20_nominal,
    Master_three_sets = dmp3
  ),
  file = "R2_1_DMP_background_comparison_summary.xlsx",
  overwrite = TRUE
)

############################################################
## 13. Final message
############################################################

message("R2-1 analysis completed.")
message("Key files generated:")
message("1. R2_1_DMP_set_counts.csv")
message("2. R2_1_DMP_direction_summary.csv")
message("3. R2_1_gene_feature_distribution.csv")
message("4. R2_1_CpG_feature_distribution.csv")
message("5. R2_1_feature_chisq_tests.txt")
message("6. R2_1_GO_enrichment_counts.csv")
message("7. R2_1_GO_top20_FDR005_all_sets.csv")
message("8. R2_1_GO_top20_nominal_all_sets.csv")
message("9. R2_1_DMP_background_comparison_summary.xlsx")


############################################################
############################################################

if (!exists("set_counts") && file.exists("R2_1_DMP_set_counts.csv")) {
  set_counts <- read.csv("R2_1_DMP_set_counts.csv")
}

if (!exists("direction_summary") && file.exists("R2_1_DMP_direction_summary.csv")) {
  direction_summary <- read.csv("R2_1_DMP_direction_summary.csv")
}

if (!exists("gene_feature_dist") && file.exists("R2_1_gene_feature_distribution.csv")) {
  gene_feature_dist <- read.csv("R2_1_gene_feature_distribution.csv")
}

if (!exists("cpg_feature_dist") && file.exists("R2_1_CpG_feature_distribution.csv")) {
  cpg_feature_dist <- read.csv("R2_1_CpG_feature_distribution.csv")
}

if (!exists("go_counts") && file.exists("R2_1_GO_enrichment_counts.csv")) {
  go_counts <- read.csv("R2_1_GO_enrichment_counts.csv")
}

if (!exists("go_top20_nominal") && file.exists("R2_1_GO_top20_nominal_all_sets.csv")) {
  go_top20_nominal <- read.csv("R2_1_GO_top20_nominal_all_sets.csv")
}

top_n_by_set <- function(df, n = 10) {
  if (!"Set" %in% names(df)) return(df)
  out <- do.call(
    rbind,
    lapply(split(df, df$Set), function(x) head(x, n))
  )
  rownames(out) <- NULL
  return(out)
}

out_file <- "R2_1_results_for_copy.txt"

capture.output({
  
  cat("============================================================\n")
  cat("R2-1 KEY RESULTS FOR COPY\n")
  cat("============================================================\n\n")
  
  cat("1. DMP set counts\n")
  cat("------------------------------------------------------------\n")
  print(set_counts)
  cat("\n")
  
  cat("2. Direction distribution\n")
  cat("------------------------------------------------------------\n")
  print(direction_summary)
  cat("\n")
  
  cat("3. Gene feature distribution\n")
  cat("------------------------------------------------------------\n")
  print(gene_feature_dist)
  cat("\n")
  
  cat("4. CpG feature distribution\n")
  cat("------------------------------------------------------------\n")
  print(cpg_feature_dist)
  cat("\n")
  
  cat("5. Chi-square tests for feature distributions\n")
  cat("------------------------------------------------------------\n")
  
  if (exists("gene_feature_mat")) {
    cat("\nGene feature contingency table:\n")
    print(gene_feature_mat)
  }
  
  if (exists("gene_feature_chisq")) {
    cat("\nGene feature chi-square test:\n")
    print(gene_feature_chisq)
  } else if (file.exists("R2_1_feature_chisq_tests.txt")) {
    cat("\nGene/CpG feature chi-square tests from txt file:\n")
    cat(paste(readLines("R2_1_feature_chisq_tests.txt"), collapse = "\n"))
    cat("\n")
  }
  
  if (exists("cpg_feature_mat")) {
    cat("\nCpG feature contingency table:\n")
    print(cpg_feature_mat)
  }
  
  if (exists("cpg_feature_chisq")) {
    cat("\nCpG feature chi-square test:\n")
    print(cpg_feature_chisq)
  }
  
  cat("\n")
  
  cat("6. GO enrichment counts by missMethyl::gometh\n")
  cat("------------------------------------------------------------\n")
  print(go_counts)
  cat("\n")
  
  cat("7. Top nominal GO terms by set, for reference only\n")
  cat("   Note: These are not FDR-significant unless FDR < 0.05.\n")
  cat("------------------------------------------------------------\n")
  
  if (exists("go_top20_nominal")) {
    cols_to_show <- intersect(
      c("Set", "ONTOLOGY", "TERM", "N", "DE", "P.DE", "FDR"),
      names(go_top20_nominal)
    )
    print(top_n_by_set(go_top20_nominal[, cols_to_show, drop = FALSE], n = 10))
  } else {
    cat("go_top20_nominal object/file not found.\n")
  }
  
  cat("\n")
  
  cat("8. Top Fisher enrichment results for gene features\n")
  cat("------------------------------------------------------------\n")
  if (exists("gene_feature_fisher")) {
    print(head(gene_feature_fisher, 20))
  } else if (file.exists("R2_1_gene_feature_fisher_enrichment.csv")) {
    tmp <- read.csv("R2_1_gene_feature_fisher_enrichment.csv")
    print(head(tmp, 20))
  } else {
    cat("Gene feature Fisher results not found.\n")
  }
  
  cat("\n")
  
  cat("9. Top Fisher enrichment results for CpG features\n")
  cat("------------------------------------------------------------\n")
  if (exists("cpg_feature_fisher")) {
    print(head(cpg_feature_fisher, 20))
  } else if (file.exists("R2_1_CpG_feature_fisher_enrichment.csv")) {
    tmp <- read.csv("R2_1_CpG_feature_fisher_enrichment.csv")
    print(head(tmp, 20))
  } else {
    cat("CpG feature Fisher results not found.\n")
  }
  
  cat("\n============================================================\n")
  cat("END OF R2-1 KEY RESULTS\n")
  cat("============================================================\n")
  
}, file = out_file)

cat(paste(readLines(out_file), collapse = "\n"))

cat("\n\nResults were also saved to: ", out_file, "\n")


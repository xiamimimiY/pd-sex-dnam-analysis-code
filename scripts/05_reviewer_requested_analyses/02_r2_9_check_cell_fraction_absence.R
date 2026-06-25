
##2026-6-17
############################################################
## Step 1: Load BL Rdata and inspect objects
############################################################

rm(list = ls())

file_bl <- "ppmi140_champ_myDMP1_BL_PD213.Rdata"

env <- new.env()
load(file_bl, envir = env)

obj_names <- ls(env)
print(obj_names)

obj_summary <- data.frame(
  object = obj_names,
  class = sapply(obj_names, function(x) paste(class(env[[x]]), collapse = "; ")),
  dim = sapply(obj_names, function(x) {
    d <- dim(env[[x]])
    if (is.null(d)) return(NA)
    paste(d, collapse = " x ")
  }),
  stringsAsFactors = FALSE
)

print(obj_summary)
write.csv(obj_summary, "BL_Rdata_object_summary.csv", row.names = FALSE)


############################################################
## Step 2: Search for estimated blood cell proportion columns
############################################################

cell_keywords <- c(
  "CD8T", "CD4T", "NK", "Bcell", "B_cell",
  "Mono", "Gran", "Neu", "Neutrophil", "Neutrophils",
  "Eos", "Bas"
)

find_cell_columns <- function(x, object_name = "object") {
  res <- list()
  
  if (is.data.frame(x) || is.matrix(x)) {
    cn <- colnames(x)
    rn <- rownames(x)
    
    col_hits <- intersect(cell_keywords, cn)
    row_hits <- intersect(cell_keywords, rn)
    
    if (length(col_hits) > 0 || length(row_hits) > 0) {
      res[[length(res) + 1]] <- data.frame(
        object = object_name,
        class = paste(class(x), collapse = "; "),
        dim = paste(dim(x), collapse = " x "),
        matched_columns = paste(col_hits, collapse = ", "),
        matched_rows = paste(row_hits, collapse = ", "),
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (is.list(x) && !is.data.frame(x)) {
    nms <- names(x)
    if (!is.null(nms)) {
      for (nm in nms) {
        res <- c(res, find_cell_columns(x[[nm]], paste0(object_name, "$", nm)))
      }
    }
  }
  
  return(res)
}

hits <- list()

for (nm in obj_names) {
  hits <- c(hits, find_cell_columns(env[[nm]], nm))
}

if (length(hits) > 0) {
  cell_hits <- do.call(rbind, hits)
} else {
  cell_hits <- data.frame()
}

print(cell_hits)
write.csv(cell_hits, "BL_Rdata_candidate_cell_objects.csv", row.names = FALSE)


############################################################
## Step 3: Check potentially useful large objects
############################################################

large_objects <- obj_summary[
  grepl("matrix|data.frame|list", obj_summary$class) & !is.na(obj_summary$dim),
]

print(large_objects)

############################################################
## Inspect myDMP1_BL_PD213 internal structure
############################################################

x <- env$myDMP1_BL_PD213

cat("Class of myDMP1_BL_PD213:\n")
print(class(x))

cat("\nLength of myDMP1_BL_PD213:\n")
print(length(x))

cat("\nNames of myDMP1_BL_PD213:\n")
print(names(x))

cat("\nStructure, max.level = 2:\n")
str(x, max.level = 2)

############################################################
## Summarize each element inside the list
############################################################

list_summary <- data.frame(
  element = if (is.null(names(x))) paste0("element_", seq_along(x)) else names(x),
  class = sapply(x, function(z) paste(class(z), collapse = "; ")),
  dim = sapply(x, function(z) {
    d <- dim(z)
    if (is.null(d)) return(NA)
    paste(d, collapse = " x ")
  }),
  length = sapply(x, length),
  stringsAsFactors = FALSE
)

print(list_summary)
write.csv(list_summary, "myDMP1_BL_PD213_list_summary.csv", row.names = FALSE)

############################################################
## Print column names of any data.frame/matrix elements
############################################################

for (i in seq_along(x)) {
  nm <- if (is.null(names(x))) paste0("element_", i) else names(x)[i]
  
  if (is.data.frame(x[[i]]) || is.matrix(x[[i]])) {
    cat("\n============================================================\n")
    cat("Element:", nm, "\n")
    cat("Class:", paste(class(x[[i]]), collapse = "; "), "\n")
    cat("Dim:", paste(dim(x[[i]]), collapse = " x "), "\n")
    cat("Column names:\n")
    print(colnames(x[[i]]))
    cat("First rows:\n")
    print(head(x[[i]], 3))
  }
}

############################################################
# RA AUTS2 preranked GSEA using the complete ranked gene list
#
# Main features:
# 1. Use RA/Treat samples only.
# 2. Split samples into AUTS2-high and AUTS2-low groups by median.
# 3. Use limma moderated t statistics for ranking.
# 4. Do not filter genes by nominal P value before GSEA.
# 5. Use KEGG gene sets in official gene-symbol format.
# 6. Use BH-adjusted P value < 0.05 as the enrichment threshold.
############################################################

rm(list = ls())
gc()

############################################################
# 1. Load packages
############################################################

required_packages <- c(
  "limma",
  "clusterProfiler",
  "enrichplot",
  "ggplot2"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Package not installed: ", pkg,
        ". Please install it before running this script."
      )
    )
  }
}

library(limma)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

############################################################
# 2. Set paths and reproducibility options
############################################################

set.seed(1234)

work_dir <- "D:/mechanical-related/GSEA/RA1"
setwd(work_dir)

expression_file <- file.path(work_dir, "merge.normalize.txt")
target_gene_file <- file.path(work_dir, "interGenes.List.txt")
gmt_file <- file.path(work_dir, "c2.cp.kegg.Hs.symbols.gmt")

output_dir <- file.path(work_dir, "AUTS2_GSEA_output")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

############################################################
# 3. Check input files
############################################################

input_files <- c(expression_file, target_gene_file, gmt_file)
missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "The following input files were not found:\n",
      paste(missing_files, collapse = "\n")
    )
  )
}

############################################################
# 4. Read and prepare expression matrix
############################################################

raw_data <- read.table(
  expression_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = ""
)

if (ncol(raw_data) < 3) {
  stop("The expression matrix contains fewer than two samples.")
}

gene_column <- colnames(raw_data)[1]
gene_names <- toupper(trimws(as.character(raw_data[[gene_column]])))

expression_data <- raw_data[, -1, drop = FALSE]
expression_data <- as.matrix(expression_data)
mode(expression_data) <- "numeric"

rownames(expression_data) <- gene_names

# Average duplicate gene symbols if any remain.
expression_data <- limma::avereps(expression_data)

# Remove genes with no positive average expression.
expression_data <- expression_data[
  rowMeans(expression_data, na.rm = TRUE) > 0,
  ,
  drop = FALSE
]

if (nrow(expression_data) == 0) {
  stop("No genes remained after expression filtering.")
}

############################################################
# 5. Keep RA disease samples and define AUTS2 groups
############################################################

sample_names <- colnames(expression_data)

# The supplied matrix uses sample names ending in _Control or _Treat.
sample_type <- sub(
  ".*_(Control|Treat)$",
  "\\1",
  sample_names
)

if (!all(sample_type %in% c("Control", "Treat"))) {
  warning(
    "Some sample names were not recognized as Control or Treat. " ,
    "Please inspect sample_type before proceeding."
  )
}

ra_data <- expression_data[, sample_type == "Treat", drop = FALSE]

if (ncol(ra_data) < 4) {
  stop("Fewer than four RA/Treat samples were detected.")
}

target_genes <- readLines(target_gene_file, warn = FALSE)
target_genes <- unique(toupper(trimws(target_genes)))
target_genes <- target_genes[target_genes != ""]

if (length(target_genes) == 0) {
  stop("No target gene was found in interGenes.List.txt.")
}

target_gene <- target_genes[1]

if (!target_gene %in% rownames(ra_data)) {
  stop(
    paste0(
      "Target gene ", target_gene,
      " was not found in the expression matrix."
    )
  )
}

target_expression <- as.numeric(ra_data[target_gene, ])
target_group <- ifelse(
  target_expression >= median(target_expression, na.rm = TRUE),
  "AUTS2_high",
  "AUTS2_low"
)

target_group <- factor(
  target_group,
  levels = c("AUTS2_low", "AUTS2_high")
)

group_counts <- table(target_group)

if (any(group_counts < 2)) {
  stop("One AUTS2 group contains fewer than two samples.")
}

write.table(
  data.frame(
    Sample = colnames(ra_data),
    AUTS2_expression = target_expression,
    AUTS2_group = target_group,
    stringsAsFactors = FALSE
  ),
  file = file.path(output_dir, "AUTS2_sample_groups.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# 6. Differential analysis between AUTS2-high and AUTS2-low
#
# All genes with valid statistics are retained for GSEA ranking.
# No nominal P-value filtering is performed at this stage.
############################################################

design <- model.matrix(~ 0 + target_group)
colnames(design) <- levels(target_group)

fit <- lmFit(ra_data, design)

contrast_matrix <- makeContrasts(
  AUTS2_high - AUTS2_low,
  levels = design
)

fit_contrast <- contrasts.fit(fit, contrast_matrix)
fit_contrast <- eBayes(fit_contrast)

limma_table <- topTable(
  fit_contrast,
  number = Inf,
  sort.by = "none",
  adjust.method = "BH"
)

limma_table$Gene <- rownames(limma_table)
limma_table <- limma_table[
  is.finite(limma_table$t) &
    !is.na(limma_table$Gene) &
    limma_table$Gene != "",
  ,
  drop = FALSE
]

limma_table <- limma_table[
  !duplicated(limma_table$Gene),
  ,
  drop = FALSE
]

limma_table <- limma_table[
  , c(
    "Gene", "logFC", "AveExpr", "t", "P.Value",
    "adj.P.Val", "B"
  ),
  drop = FALSE
]

write.table(
  limma_table,
  file = file.path(output_dir, "AUTS2_high_vs_low_limma_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# 7. Build the complete ranked gene list
############################################################

# Ranking statistic: limma moderated t statistic.
# All eligible genes are used; no P-value filtering is applied.
gene_list <- limma_table$t
names(gene_list) <- limma_table$Gene
gene_list <- gene_list[is.finite(gene_list)]
gene_list <- sort(gene_list, decreasing = TRUE)

write.table(
  data.frame(
    Gene = names(gene_list),
    Ranking_statistic = as.numeric(gene_list),
    stringsAsFactors = FALSE
  ),
  file = file.path(output_dir, "AUTS2_complete_ranked_gene_list.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# 8. Read KEGG gene sets and run preranked GSEA
############################################################

kegg_sets <- read.gmt(gmt_file)

gsea_object <- GSEA(
  geneList = gene_list,
  TERM2GENE = kegg_sets,
  exponent = 1,
  minGSSize = 5,
  maxGSSize = 5000,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  seed = TRUE,
  verbose = FALSE
)

gsea_table <- as.data.frame(gsea_object)

if (nrow(gsea_table) == 0) {
  stop("No gene sets were returned by GSEA.")
}

# Complete enrichment table: no pathway filtering.
write.table(
  gsea_table,
  file = file.path(output_dir, "AUTS2_GSEA_full_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Significant enrichment table: BH-FDR < 0.05.
gsea_fdr <- gsea_table[
  !is.na(gsea_table$p.adjust) &
    gsea_table$p.adjust < 0.05,
  ,
  drop = FALSE
]

write.table(
  gsea_fdr,
  file = file.path(output_dir, "AUTS2_GSEA_FDR05_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# 9. Plot up to five FDR-significant pathways in each direction
############################################################

select_terms <- function(tab, direction, n = 5) {
  selected <- tab[
    if (direction == "high") tab$NES > 0 else tab$NES < 0,
    ,
    drop = FALSE
  ]

  selected <- selected[order(selected$p.adjust, -abs(selected$NES)), , drop = FALSE]
  head(selected$ID, n)
}

high_terms <- select_terms(gsea_fdr, "high", 5)
low_terms <- select_terms(gsea_fdr, "low", 5)

if (length(high_terms) > 0) {
  pdf(
    file.path(output_dir, "AUTS2_high_enriched_GSEA.pdf"),
    width = 7,
    height = 8
  )
  print(
    gseaplot2(
      gsea_object,
      geneSetID = high_terms,
      title = "Pathways enriched in AUTS2-high RA samples",
      base_size = 9
    )
  )
  dev.off()
}

if (length(low_terms) > 0) {
  pdf(
    file.path(output_dir, "AUTS2_low_enriched_GSEA.pdf"),
    width = 7,
    height = 8
  )
  print(
    gseaplot2(
      gsea_object,
      geneSetID = low_terms,
      title = "Pathways enriched in AUTS2-low RA samples",
      base_size = 9
    )
  )
  dev.off()
}

############################################################
# 10. Save a compact analysis summary
############################################################

summary_lines <- c(
  paste0("Target gene: ", target_gene),
  paste0("RA/Treat samples: ", ncol(ra_data)),
  paste0("AUTS2-low samples: ", unname(group_counts["AUTS2_low"])),
  paste0("AUTS2-high samples: ", unname(group_counts["AUTS2_high"])),
  paste0("Genes in complete ranked list: ", length(gene_list)),
  paste0("All tested gene sets: ", nrow(gsea_table)),
  paste0("FDR < 0.05 gene sets: ", nrow(gsea_fdr)),
  "Ranking statistic: limma moderated t statistic",
  "Pre-ranking P-value filter: none",
  "Gene identifiers: official gene symbols",
  "Gene-set file: c2.cp.kegg.Hs.symbols.gmt",
  "Significance threshold: Benjamini-Hochberg adjusted P value < 0.05",
  "Random seed: 1234"
)

writeLines(
  summary_lines,
  con = file.path(output_dir, "AUTS2_GSEA_analysis_summary.txt")
)

cat("GSEA completed successfully.\n")
cat("Results saved to:\n", output_dir, "\n")

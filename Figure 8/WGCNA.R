#!/usr/bin/env Rscript

# T cell hdWGCNA workflow.
# Input defaults to ../7pbmc_celltype.rds when this script is run from results3.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(WGCNA)
  library(hdWGCNA)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# ---- 1. Paths and parameters ----
seurat_obj_path <- "D:/mechanical-related/7-SLE-scRNA/24/WGCNA/7pbmc_celltype.rds"
if (!file.exists(seurat_obj_path) && file.exists("7pbmc_celltype.rds")) {
  seurat_obj_path <- "7pbmc_celltype.rds"
}

outdir <- file.path(getwd(), "Tcell_WGCNA")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

celltype_col <- "celltype.main"
target_celltype <- "T_cell"
donor_col <- "sample"
assay_use <- "RNA"
layer_use <- "data"
wgcna_name <- "Tcell"
k <- 20
max_shared <- 10
n_threads <- 10
gene_select_fraction <- 0.05
marker_gene <- "AUTS2"

# ---- 2. Helpers ----
get_assay_data_compat <- function(obj, assay = "RNA", layer = "data") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
}

save_plot_pdf <- function(path, plot, width = 8, height = 5) {
  ggsave(path, plot = plot, width = width, height = height)
}

# ---- 3. Read data and subset T cells ----
if (!file.exists(seurat_obj_path)) {
  stop("Cannot find Seurat object: ", seurat_obj_path, call. = FALSE)
}

cat(">>> Reading Seurat object:", seurat_obj_path, "\n")
seurat_obj <- readRDS(seurat_obj_path)

if (!celltype_col %in% colnames(seurat_obj@meta.data)) {
  stop("Metadata column not found: ", celltype_col, call. = FALSE)
}

cat(">>> Metadata columns:\n")
print(colnames(seurat_obj@meta.data))

if (!donor_col %in% colnames(seurat_obj@meta.data)) {
  donor_candidates <- unique(c(
    donor_col,
    "sample", "Sample", "sample_id", "SampleID", "orig.ident",
    "donor", "Donor", "patient", "Patient", "subject", "Subject",
    "individual", "Individual", "case", "Case"
  ))
  donor_hit <- donor_candidates[donor_candidates %in% colnames(seurat_obj@meta.data)]

  if (length(donor_hit) > 0) {
    donor_col <- donor_hit[[1]]
    cat(">>> Metadata column 'sample' was not found; using donor column:", donor_col, "\n")
  } else {
    donor_col <- "donor"
    seurat_obj$donor <- "all_cells"
    cat(">>> No donor/sample metadata column found; created donor = 'all_cells'.\n")
  }
}

cat(">>> Cell types in", celltype_col, ":\n")
print(table(seurat_obj@meta.data[[celltype_col]]))

if (!target_celltype %in% seurat_obj@meta.data[[celltype_col]]) {
  stop("Target cell type not found: ", target_celltype, call. = FALSE)
}

cat(">>> Subsetting T cells...\n")
tcell_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[celltype_col]] == target_celltype]
tcell_obj <- subset(seurat_obj, cells = tcell_cells)
tcell_obj$cell_type <- target_celltype
tcell_obj$donor <- tcell_obj@meta.data[[donor_col]]

cat(">>> T cell number:", ncol(tcell_obj), "\n")
cat(">>> Donor distribution:\n")
print(table(tcell_obj$donor))

saveRDS(tcell_obj, file.path(outdir, "Tcell_subset.rds"))

# ---- 4. Build metacells ----
cat(">>> Setting up hdWGCNA object...\n")
tcell_obj <- SetupForWGCNA(
  tcell_obj,
  gene_select = "fraction",
  fraction = gene_select_fraction,
  wgcna_name = wgcna_name
)

cat(">>> Constructing metacells, k =", k, "...\n")
tcell_obj <- MetacellsByGroups(
  seurat_obj = tcell_obj,
  group.by = c("cell_type", "donor"),
  k = k,
  max_shared = max_shared,
  ident.group = "cell_type"
)

tcell_obj <- NormalizeMetacells(tcell_obj)
metacell_obj <- GetMetacellObject(tcell_obj)

cat(">>> Metacell object dimensions:\n")
print(dim(metacell_obj))
cat(">>> Metacell donor distribution:\n")
print(table(metacell_obj$donor))

saveRDS(tcell_obj, file.path(outdir, "Tcell_metacell.rds"))

# ---- 5. Set expression matrix for WGCNA ----
cat(">>> Setting expression data for WGCNA...\n")
tcell_obj <- SetDatExpr(
  tcell_obj,
  group_name = target_celltype,
  group.by = "cell_type",
  assay = assay_use,
  slot = layer_use
)

# ---- 6. Soft-power selection ----
cat(">>> Testing soft powers...\n")
tcell_obj <- TestSoftPowers(tcell_obj, networkType = "signed")

plot_list <- PlotSoftPowers(tcell_obj)
p_soft <- wrap_plots(plot_list, ncol = 2)
save_plot_pdf(file.path(outdir, "01_soft_power.pdf"), p_soft, width = 10, height = 6)

power_table <- GetPowerTable(tcell_obj)
write_csv(power_table, file.path(outdir, "soft_power_table.csv"))

candidate_power <- power_table %>%
  filter(SFT.R.sq >= 0.8) %>%
  pull(Power)

if (length(candidate_power) > 0) {
  soft_power <- min(candidate_power, na.rm = TRUE)
} else {
  soft_power <- power_table %>%
    arrange(desc(SFT.R.sq)) %>%
    slice(1) %>%
    pull(Power)
  cat(">>> No power reached SFT.R.sq >= 0.8; using best available power.\n")
}

cat(">>> Selected soft_power =", soft_power, "\n")

# ---- 7. Construct co-expression network ----
cat(">>> Constructing T cell co-expression network...\n")
tcell_obj <- ConstructNetwork(
  tcell_obj,
  soft_power = soft_power,
  setDatExpr = FALSE,
  tom_name = wgcna_name,
  tom_outdir = outdir,
  overwrite_tom = TRUE,
  nThreads = n_threads
)

pdf(file.path(outdir, "02_dendrogram.pdf"), width = 8, height = 5)
PlotDendrogram(tcell_obj, main = "T cell hdWGCNA Dendrogram")
dev.off()

# ---- 8. Module eigengenes and connectivity ----
cat(">>> Calculating module eigengenes...\n")
tcell_obj <- ModuleEigengenes(
  tcell_obj,
  group.by.vars = "donor"
)

hMEs <- GetMEs(tcell_obj, harmonized = TRUE)
if ("grey" %in% colnames(hMEs)) {
  hMEs <- hMEs[, colnames(hMEs) != "grey", drop = FALSE]
}
write_csv(
  as.data.frame(hMEs) %>% rownames_to_column("cell"),
  file.path(outdir, "harmonized_module_eigengenes.csv")
)

cat(">>> Calculating module connectivity...\n")
tcell_obj <- ModuleConnectivity(
  tcell_obj,
  group.by = "cell_type",
  group_name = target_celltype
)

p_kme <- PlotKMEs(tcell_obj, ncol = 5)
save_plot_pdf(file.path(outdir, "03_KME_plots.pdf"), p_kme, width = 12, height = 6)

modules <- GetModules(tcell_obj)
write_csv(modules, file.path(outdir, "modules.csv"))
cat(">>> Module composition:\n")
print(table(modules$module))

# ---- 9. Optional marker-module correlation ----
if (marker_gene %in% rownames(tcell_obj) && marker_gene %in% rownames(metacell_obj) && ncol(hMEs) > 0) {
  cat(">>> Calculating", marker_gene, "correlation with module eigengenes...\n")

  marker_expr <- get_assay_data_compat(metacell_obj, assay = assay_use, layer = layer_use)[marker_gene, ] %>%
    as.numeric()
  marker_df <- tibble(cell = colnames(metacell_obj), marker_expr = marker_expr)

  me_df <- hMEs %>%
    as.data.frame() %>%
    rownames_to_column("cell")

  cor_df <- marker_df %>% left_join(me_df, by = "cell")

  cor_res <- map_dfr(colnames(hMEs), function(mod) {
    x <- cor_df$marker_expr
    y <- cor_df[[mod]]
    valid <- complete.cases(x, y)
    if (sum(valid) < 5 || sd(x[valid]) == 0 || sd(y[valid]) == 0) {
      return(tibble(module = mod, estimate = NA_real_, p_value = NA_real_))
    }
    test <- cor.test(x[valid], y[valid], method = "pearson")
    tibble(module = mod, estimate = unname(test$estimate), p_value = test$p.value)
  }) %>%
    filter(!is.na(p_value)) %>%
    mutate(
      fdr = p.adjust(p_value, method = "fdr"),
      signif = case_when(
        fdr < 0.001 ~ "***",
        fdr < 0.01 ~ "**",
        fdr < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )

  write_csv(cor_res, file.path(outdir, paste0(marker_gene, "_module_correlation.csv")))

  p_bar <- cor_res %>%
    mutate(module = forcats::fct_reorder(module, estimate)) %>%
    ggplot(aes(x = module, y = estimate, fill = estimate)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = signif), vjust = ifelse(cor_res$estimate > 0, -0.3, 1.3), size = 3) +
    scale_fill_gradient2(low = "#2187A5", mid = "white", high = "#EF4028", midpoint = 0) +
    labs(
      title = paste(marker_gene, "correlation with T cell module eigengenes"),
      x = "Module",
      y = "Pearson r"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))

  save_plot_pdf(
    file.path(outdir, paste0(marker_gene, "_module_correlation_barplot.pdf")),
    p_bar,
    width = 8,
    height = 5
  )
} else {
  cat(">>> Skip marker-module correlation: marker missing or no non-grey modules.\n")
}

# ---- 10. Save final object ----
saveRDS(tcell_obj, file.path(outdir, "Tcell_wgcna.rds"))

cat("\n>>> ====== T cell WGCNA finished ======\n")
cat(">>> Results saved in:", outdir, "\n")
cat(">>> Key outputs:\n")
cat("    - Tcell_subset.rds\n")
cat("    - Tcell_metacell.rds\n")
cat("    - Tcell_wgcna.rds\n")
cat("    - modules.csv\n")
cat("    - soft_power_table.csv\n")
cat("    - 01_soft_power.pdf\n")
cat("    - 02_dendrogram.pdf\n")
cat("    - 03_KME_plots.pdf\n")

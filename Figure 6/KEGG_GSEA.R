# Run with Source in RStudio or with Rscript.
# Reuses the project's donor-paired Hallmark ranking method with KEGG Legacy sets.
local({
  analysis_dir <- "D:/mechanical-related/111"
  seurat_file <- file.path(analysis_dir, "Fibroblast_celltype.rds")
  kegg_gmt <- file.path(dirname(analysis_dir),
                       "GSEA", "RA", "c2.cp.kegg.Hs.symbols.gmt")
  output_dir <- file.path(analysis_dir, "AUTS2_high_low_KEGG_results")
  fdr_cutoff <- 0.05
  max_pathways_per_direction <- 10L
  min_cells_per_group <- 20L

  packages <- c("Seurat", "Matrix", "edgeR", "limma", "fgsea", "BiocParallel", "ggplot2")
  available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  if (any(!available)) stop("Missing R packages: ", paste(packages[!available], collapse = ", "))
  for (path in c(seurat_file, kegg_gmt)) {
    if (!file.exists(path)) stop("Input file not found: ", path)
  }

  message("Reading fibroblast RNA counts and assigning AUTS2 groups...")
  cells <- readRDS(seurat_file)
  counts <- Seurat::GetAssayData(cells, assay = "RNA", layer = "counts")
  normalized <- Seurat::GetAssayData(cells, assay = "RNA", layer = "data")
  if (!"AUTS2" %in% rownames(normalized) || !"orig.ident" %in% names(cells[[]])) {
    stop("Expected the AUTS2 gene in RNA data and donor IDs in orig.ident.")
  }
  expression <- as.numeric(normalized["AUTS2", colnames(counts)])
  cutoff <- median(expression)
  group_levels <- c("AUTS2_low", "AUTS2_high")
  group <- factor(ifelse(expression > cutoff, group_levels[2L], group_levels[1L]),
                  levels = group_levels)
  donor <- as.character(cells[[]][colnames(counts), "orig.ident"])
  if (anyNA(donor) || any(!nzchar(donor))) stop("Some cells have no donor ID.")
  group_counts <- table(donor, group)
  eligible_donors <- rownames(group_counts)[apply(group_counts, 1L, min) >= min_cells_per_group]
  if (length(eligible_donors) < 3L) stop("Fewer than three donors have enough cells in both groups.")
  donor_summary <- data.frame(donor = rownames(group_counts),
                              AUTS2_low = as.integer(group_counts[, "AUTS2_low"]),
                              AUTS2_high = as.integer(group_counts[, "AUTS2_high"]))
  donor_summary$included <- donor_summary$donor %in% eligible_donors
  message("AUTS2 median: ", cutoff, "; high: ", sum(group == "AUTS2_high"),
          "; low: ", sum(group == "AUTS2_low"), "; paired donors: ", length(eligible_donors))

  # One pseudobulk library per donor and group, with the donor as a paired block.
  samples <- expand.grid(group = group_levels, donor = eligible_donors,
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  rownames(samples) <- paste(samples$donor, samples$group, sep = "__")
  sample_index <- match(paste(donor, group, sep = "__"), rownames(samples))
  included_cells <- which(!is.na(sample_index))
  membership <- Matrix::sparseMatrix(i = included_cells, j = sample_index[included_cells],
                                     x = 1, dims = c(ncol(counts), nrow(samples)))
  pseudobulk <- as.matrix(counts %*% membership)
  colnames(pseudobulk) <- rownames(samples)
  rm(cells, counts, normalized, membership)
  samples$group <- factor(samples$group, levels = group_levels)
  samples$donor <- factor(samples$donor, levels = eligible_donors)
  design <- stats::model.matrix(~ donor + group, data = samples)

  # Omit the grouping gene and rank all remaining expressed genes by moderated t.
  message("Fitting the paired donor model...")
  y <- edgeR::DGEList(counts = pseudobulk[rownames(pseudobulk) != "AUTS2", , drop = FALSE])
  expressed <- edgeR::filterByExpr(y, design = design)
  y <- edgeR::calcNormFactors(y[expressed, , keep.lib.sizes = FALSE], method = "TMM")
  fit <- limma::eBayes(limma::lmFit(limma::voom(y, design, plot = FALSE), design))
  t_stat <- fit$t[, "groupAUTS2_high"]
  t_stat <- t_stat[is.finite(t_stat)]
  ranks <- t_stat[order(-t_stat, names(t_stat), method = "radix")]

  gmt_rows <- strsplit(readLines(kegg_gmt, warn = FALSE), "\t", fixed = TRUE)
  gmt_rows <- gmt_rows[lengths(gmt_rows) >= 3L]
  pathways <- setNames(lapply(gmt_rows, function(row) setdiff(row[-c(1L, 2L)], "AUTS2")),
                       vapply(gmt_rows, function(row) row[[1L]], character(1)))
  if (!length(pathways) || anyDuplicated(names(pathways))) stop("Invalid KEGG pathway GMT.")

  message("Running KEGG Legacy GSEA with ", length(ranks), " ranked genes...")
  set.seed(1234L)
  enrichment <- as.data.frame(fgsea::fgseaMultilevel(
    pathways = pathways, stats = ranks, minSize = 10L, maxSize = 500L,
    eps = 0, nPermSimple = 10000L, BPPARAM = BiocParallel::SerialParam()
  ))
  # Adjust across all tested pathways, including any for which NES or P is unavailable.
  enrichment$padj <- stats::p.adjust(enrichment$pval, method = "BH", n = nrow(enrichment))
  enrichment <- enrichment[order(enrichment$padj, -abs(enrichment$NES), enrichment$pathway), ]
  enrichment$leadingEdge <- vapply(enrichment$leadingEdge, paste, character(1), collapse = ";")
  enrichment$direction <- ifelse(enrichment$NES > 0, "AUTS2_high",
                                 ifelse(enrichment$NES < 0, "AUTS2_low", NA_character_))
  significant <- subset(enrichment, is.finite(NES) & is.finite(padj) & padj < fdr_cutoff)
  shown <- rbind(head(subset(significant, NES > 0), max_pathways_per_direction),
                 head(subset(significant, NES < 0), max_pathways_per_direction))
  shown <- shown[order(shown$NES), ]

  title <- "AUTS2_high vs AUTS2_low"
  subtitle <- sprintf("KEGG Legacy GSEA | %d paired donors | BH FDR < %g",
                      length(eligible_donors), fdr_cutoff)
  if (nrow(shown) > 0L) {
    pathway_labels <- gsub("_", " ", sub("^KEGG_", "", shown$pathway))
    pathway_labels <- vapply(pathway_labels, function(label) {
      paste(strwrap(label, width = 45L), collapse = "\n")
    }, character(1))
    shown$label <- factor(pathway_labels, levels = unique(pathway_labels))
    plot <- ggplot2::ggplot(shown, ggplot2::aes(x = NES, y = label, fill = direction)) +
      ggplot2::geom_vline(xintercept = 0, color = "grey60", linewidth = 0.35) +
      ggplot2::geom_col(width = 0.72) +
      ggplot2::scale_fill_manual(values = c(AUTS2_low = "#1C7C99", AUTS2_high = "#E84426"),
                                 breaks = group_levels, name = "Enriched in") +
      ggplot2::labs(title = title, subtitle = subtitle, x = "Normalized enrichment score (NES)",
                    y = NULL) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5))
  } else {
    plot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No KEGG pathways passed the FDR cutoff.") +
      ggplot2::labs(title = title, subtitle = subtitle) +
      ggplot2::theme_void(base_size = 12)
  }

  if (!dir.exists(output_dir)) dir.create(output_dir)
  if (!dir.exists(output_dir)) stop("Could not create output folder: ", output_dir)
  write.csv(donor_summary, file.path(output_dir, "AUTS2_groups_by_donor.csv"), row.names = FALSE)
  write.csv(enrichment, file.path(output_dir, "KEGG_GSEA_all_pathways.csv"), row.names = FALSE)
  pdf_file <- file.path(output_dir, "AUTS2_high_low_KEGG_GSEA.pdf")
  ggplot2::ggsave(pdf_file, plot, width = 10, height = max(3.5, 1.8 + 0.34 * nrow(shown)))
  message("KEGG pathways at BH FDR < ", fdr_cutoff, ": ", nrow(significant),
          " (high: ", sum(significant$NES > 0), "; low: ", sum(significant$NES < 0), ").")
  message("Saved: ", pdf_file)
})

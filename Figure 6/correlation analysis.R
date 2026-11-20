# Run with Source in RStudio or with Rscript.
local({
  analysis_dir <- "D:/mechanical-related/111"
  output_dir <- file.path(analysis_dir, "RA_AUTS2_Mechano_BH_reanalysis")
  object_file <- file.path(analysis_dir, "Fibroblast_celltype.rds")
  # Keep the same module order as the RA correlation reanalysis for reproducible scoring.
  modules <- c(
    "Actin cytoskeleton", "Actin filament based movement", "Actin filament polymerization",
    "Cell-cell adhesion", "Cellular response to mechanical stimulus", "Extracellular matrix assembly",
    "Extracellular matrix constituent secretion", "Extracellular matrix disassembly",
    "Extracellular matrix organization", "Hippo signaling", "Integrin activation", "Microtubule",
    "Protein complex involved in cell adhesion", "Response to mechanical stimulus",
    "Rho protein signal transduction"
  )
  module_ids <- gsub("[^A-Za-z0-9_]", "_", modules)
  gene_files <- file.path(analysis_dir, paste0(modules, ".txt"))
  for (path in c(object_file, gene_files)) {
    if (!file.exists(path)) stop("Input file not found: ", path)
  }
  library(Seurat)
  library(ggplot2)

  message("Reading RA fibroblast data and extracting AUTS2 RNA expression...")
  fls <- readRDS(object_file)
  rna <- Seurat::GetAssayData(fls, assay = "RNA", layer = "data")
  metadata <- fls[[]]
  if (!"AUTS2" %in% rownames(rna) || !"orig.ident" %in% names(metadata)) {
    stop("AUTS2 RNA expression and orig.ident donor IDs are required.")
  }
  cell_ids <- colnames(rna)
  expression <- as.numeric(rna["AUTS2", ])
  donor <- as.character(metadata[cell_ids, "orig.ident"])
  if (any(!is.finite(expression)) || anyNA(donor) || any(!nzchar(donor))) {
    stop("Some cells lack finite AUTS2 expression or a donor ID.")
  }
  donor_ids <- sort(unique(donor))
  if (length(donor_ids) != 4L) stop("Expected the four donors in the RA fibroblast object.")

  measured <- lapply(gene_files, function(path) {
    genes <- toupper(trimws(readLines(path, warn = FALSE, encoding = "UTF-8")))
    intersect(genes[nzchar(genes) & !startsWith(genes, "#")], rownames(rna))
  })
  features <- lapply(measured, setdiff, y = "AUTS2")
  control_pool <- setdiff(rownames(rna), "AUTS2")
  if (any(lengths(features) < 2L)) stop("At least two measured genes are required per module.")

  message("Scoring 15 gene sets with AUTS2 excluded...")
  score_prefix <- "MechWithoutAUTS2_"
  scored <- Seurat::AddModuleScore(
    fls, features = features, pool = control_pool, assay = "RNA", slot = "data",
    nbin = 24L, ctrl = 100L, seed = 1L, search = FALSE, name = score_prefix
  )
  score <- scored[[]][cell_ids, paste0(score_prefix, seq_along(modules)), drop = FALSE]
  names(score) <- module_ids
  stopifnot(all(is.finite(as.matrix(score))))

  cell_data <- data.frame(cell = cell_ids, donor = donor, AUTS2_expression = expression,
                          score, row.names = NULL)
  donor_data <- stats::aggregate(
    cell_data[, c("AUTS2_expression", module_ids)], by = list(donor = donor), FUN = mean
  )
  donor_data$n_cells <- as.integer(table(donor)[donor_data$donor])

  # Enumerate all 4! permutations. Repeated rank values retain their multiplicity.
  permutations <- function(indices) {
    if (length(indices) == 1L) return(matrix(indices, nrow = 1L))
    do.call(rbind, lapply(indices, function(i) {
      cbind(i, permutations(indices[indices != i]))
    }))
  }
  perm_indices <- permutations(seq_along(donor_ids))
  stopifnot(nrow(perm_indices) == factorial(length(donor_ids)))
  donor_test <- function(x, y) {
    if (length(unique(x)) < 2L || length(unique(y)) < 2L) return(c(rho = NA, p_val = NA))
    rho <- cor(x, y, method = "spearman")
    null <- apply(perm_indices, 1L, function(i) cor(x, y[i], method = "spearman"))
    c(rho = rho, p_val = mean(abs(null) >= abs(rho) - 1e-12))
  }
  cell_test <- function(x, y) {
    if (length(unique(x)) < 2L || length(unique(y)) < 2L) return(c(rho = NA, p_val = NA))
    result <- cor.test(x, y, method = "spearman", exact = FALSE)
    c(rho = unname(result$estimate), p_val = result$p.value)
  }
  summarize_tests <- function(data, test) {
    estimates <- t(vapply(module_ids, function(id) test(data$AUTS2_expression, data[[id]]),
                          c(rho = 0, p_val = 0)))
    result <- data.frame(
      Module = module_ids, Pathway = modules, n = nrow(data), estimates,
      measured_genes = lengths(measured), scored_genes = lengths(features),
      AUTS2_removed = vapply(measured, function(genes) "AUTS2" %in% genes, logical(1)),
      row.names = NULL
    )
    result$padj <- p.adjust(result$p_val, method = "BH", n = length(modules))
    result$neg_log10_padj <- -log10(pmax(result$padj, .Machine$double.xmin))
    result$significant_BH <- !is.na(result$padj) & result$padj < 0.05
    result
  }
  cell_results <- summarize_tests(cell_data, cell_test)
  donor_results <- summarize_tests(donor_data, donor_test)

  plot_data <- cell_results[order(cell_results$rho), ]
  plot_data$Pathway <- factor(plot_data$Pathway, levels = plot_data$Pathway)
  rho_limit <- max(abs(plot_data$rho), na.rm = TRUE)
  plot <- ggplot(plot_data, aes(x = 1, y = Pathway)) +
    geom_point(aes(size = neg_log10_padj, fill = rho), shape = 21, colour = "black",
               stroke = 0.6, alpha = 0.9, na.rm = TRUE) +
    scale_size_continuous(range = c(3, 11), name = "-log10(padj)") +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         limits = c(-rho_limit, rho_limit), name = "Spearman rho") +
    scale_y_discrete(drop = FALSE) +
    xlim(0.5, 1.5) +
    labs(title = "AUTS2 and mechanical pathways in RA fibroblasts",
         subtitle = sprintf("Exploratory cell-level analysis | %s cells from %d donors",
                            format(nrow(cell_data), big.mark = ","), nrow(donor_data)),
         caption = "AUTS2 excluded from scoring. P values adjusted across 15 pathways (BH).",
         x = NULL, y = NULL) +
    guides(fill = guide_colourbar(order = 1, barwidth = 0.6, barheight = 10),
           size = guide_legend(order = 2, override.aes = list(fill = "grey50"))) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(colour = "black"),
          panel.grid = element_blank(), legend.position = "right",
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0))

  if (!dir.exists(output_dir)) dir.create(output_dir)
  if (!dir.exists(output_dir)) stop("Could not create output folder: ", output_dir)
  write.csv(cell_results, file.path(output_dir, "RA_cell_level_correlations.csv"), row.names = FALSE)
  write.csv(donor_results, file.path(output_dir, "RA_donor_level_correlations.csv"), row.names = FALSE)
  write.csv(cell_data, file.path(output_dir, "RA_cell_scores.csv"), row.names = FALSE)
  write.csv(donor_data, file.path(output_dir, "RA_donor_mean_scores.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
  ggsave(file.path(output_dir, "RA_AUTS2_Mechano_BH_bubble.pdf"), plot, width = 6.8, height = 7)
  ggsave(file.path(output_dir, "RA_AUTS2_Mechano_BH_bubble.png"), plot,
         width = 6.8, height = 7, dpi = 600, bg = "white")
  message("Cell-level BH < 0.05: ", sum(cell_results$significant_BH), "/15")
  message("Donor-level BH < 0.05: ", sum(donor_results$significant_BH), "/15")
  message("Saved in: ", output_dir)
})

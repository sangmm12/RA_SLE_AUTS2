# Run with Source in RStudio or with Rscript.
local({
  analysis_dir <- "D:/mechanical-related/6RA/2/\u76F8\u5173\u6027"
  output_dir <- file.path(analysis_dir, "RA_correlation_padj_results")
  suppressPackageStartupMessages({
    library(linkET)
    library(ggplot2)
  })

  heat_genes <- c("CCL2", "MMP9", "TNFSF10", "ABCA1", "TOP2A", "CXCL10", "TYMS",
                  "COL17A1", "IL1RN", "OLR1", "VNN1", "MX1", "MX2", "GBP1",
                  "IFIT5", "TGFBR3", "TPM2")
  expr <- t(read.csv(file.path(analysis_dir, "Sample Type Matrix.csv"),
                     row.names = 1, check.names = FALSE))
  colnames(expr) <- trimws(colnames(expr))
  selected <- c("AUTS2", heat_genes)
  missing <- setdiff(selected, colnames(expr))
  if (length(missing)) stop("Missing genes: ", paste(missing, collapse = ", "))
  if (anyDuplicated(colnames(expr)[colnames(expr) %in% selected])) {
    stop("Duplicate gene names after trimming whitespace.")
  }
  expr <- expr[, selected, drop = FALSE]
  if (!is.numeric(expr) || any(!is.finite(expr)) || any(apply(expr, 2, sd) == 0)) {
    stop("Selected genes must have finite, numeric and varying expression.")
  }
  heat_mat <- expr[, heat_genes, drop = FALSE]

  # Keep the original sample set and linkET Mantel method. Adjust all 17 links.
  mantel <- mantel_test(spec = expr, env = heat_mat, spec_select = list(AUTS2 = 1L),
                        seed = 123L, permutations = 999L)
  stopifnot(nrow(mantel) == length(heat_genes), all(is.finite(mantel$p)))
  mantel$padj <- p.adjust(mantel$p, method = "BH")
  mantel$rd <- cut(mantel$r, c(-Inf, 0.2, 0.3, Inf), right = FALSE,
                   labels = c("< 0.2", "0.2 - < 0.3", ">= 0.3"))
  mantel$pd <- cut(mantel$padj, c(-Inf, 0.005, 0.01, 0.05, Inf), right = FALSE,
                   labels = c("< 0.005", "0.005 - < 0.01", "0.01 - < 0.05", ">= 0.05"))

  # Adjust the 136 unique Pearson pairs. The diagonal and mirrored pairs are omitted.
  pearson <- correlate(heat_mat, method = "pearson")
  upper <- upper.tri(pearson$p)
  pairs <- which(upper, arr.ind = TRUE)
  pearson_table <- data.frame(
    gene1 = rownames(pearson$r)[pairs[, "row"]],
    gene2 = colnames(pearson$r)[pairs[, "col"]],
    r = pearson$r[upper], p = pearson$p[upper]
  )
  pearson_table$padj <- p.adjust(pearson_table$p, method = "BH")
  pearson$p[upper] <- pearson_table$padj
  pearson$p[lower.tri(pearson$p)] <- t(pearson$p)[lower.tri(pearson$p)]
  diag(pearson$p) <- NA_real_

  links <- geom_couple(data = mantel, aes(colour = pd, size = rd),
                       label.size = 3.88, label.fontface = 1, nudge_x = 0.2,
                       curvature = nice_curvature(by = "from"))
  theme_and_guides <- list(
    scale_size_manual(values = setNames(c(0.5, 1.5, 3), levels(mantel$rd))),
    scale_colour_manual(values = setNames(color_pal(4, alpha = 0.6), levels(mantel$pd)),
                        drop = FALSE),
    guides(size = guide_legend(title = "Mantel's r", order = 2),
           colour = guide_legend(title = "Mantel padj", order = 1),
           fill = guide_colourbar(title = "Pearson's r", order = 3)),
    theme(legend.box.spacing = grid::unit(0, "pt"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  )
  pearson_scale <- function(palette) {
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, palette)),
                         limits = c(-1, 1), breaks = seq(-1, 1, 0.5))
  }
  square <- qcorrplot(pearson, type = "upper", diag = FALSE, grid_col = "grey50",
                      grid_size = 0.5) +
    geom_square() +
    geom_mark(size = 4, only_mark = TRUE, sig_level = c(0.05, 0.01, 0.001),
              sig_thres = 0.05, colour = "white") +
    links + pearson_scale("RdBu") + theme_and_guides +
    labs(caption = "Heatmap stars use Pearson BH padj across 136 unique gene pairs.")
  bubble <- qcorrplot(pearson, type = "upper", diag = FALSE, grid_col = NA) +
    geom_point(shape = 21, size = 8, fill = NA, stroke = 0.35, colour = "black") +
    geom_point(aes(size = abs(r), fill = r), shape = 21, stroke = 0.35, colour = "black") +
    scale_size_continuous(range = c(1, 8), guide = "none") +
    ggnewscale::new_scale("size") +
    links + pearson_scale("Spectral") + theme_and_guides

  if (!dir.exists(output_dir)) dir.create(output_dir)
  if (!dir.exists(output_dir)) stop("Could not create output folder: ", output_dir)
  write.csv(mantel, file.path(output_dir, "RA_AUTS2_Mantel_padj.csv"), row.names = FALSE)
  write.csv(pearson_table, file.path(output_dir, "RA_Pearson_gene_pairs_padj.csv"), row.names = FALSE)
  ggsave(file.path(output_dir, "RA_Gene_Corr_Mantel_Square_padj.pdf"), square, width = 8.8, height = 6)
  ggsave(file.path(output_dir, "RA_Gene_Corr_Mantel_Bubble_padj.pdf"), bubble, width = 8.8, height = 6)
  message("Samples: ", nrow(expr))
  message("Mantel P < 0.05: ", sum(mantel$p < 0.05), "; BH padj < 0.05: ", sum(mantel$padj < 0.05))
  message("Pearson P < 0.05: ", sum(pearson_table$p < 0.05),
          "; BH padj < 0.05: ", sum(pearson_table$padj < 0.05))
})

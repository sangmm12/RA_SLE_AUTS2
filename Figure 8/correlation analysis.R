# Redraw the existing cell-level Spearman results using BH-adjusted P values.
# Run with Source in RStudio or with Rscript.
local({
  library(ggplot2)

  analysis_dir <- "D:/mechanical-related/7-SLE-scRNA/24"
  input_file <- file.path(analysis_dir, "AUTS2_15pathway_Bubble_Vertical_Data.csv")
  output_dir <- file.path(analysis_dir, "AUTS2_15pathway_Bubble_padj_results")
  output_name <- "AUTS2_15pathway_Bubble_Vertical_padj"
  cor_df <- read.csv(input_file, stringsAsFactors = FALSE)
  required <- c("Module", "rho", "p_val", "ModuleSize")
  if (!all(required %in% names(cor_df)) || nrow(cor_df) != 15L ||
      anyNA(cor_df$Module) || any(!nzchar(cor_df$Module)) || anyDuplicated(cor_df$Module) ||
      any(!is.finite(cor_df$rho)) || any(!is.finite(cor_df$p_val)) ||
      any(cor_df$p_val < 0 | cor_df$p_val > 1)) {
    stop("Expected 15 unique pathways with finite Spearman rho and valid raw P values.")
  }

  # Adjust across all 15 tests before plotting any pathway.
  cor_df$padj <- p.adjust(cor_df$p_val, method = "BH")
  cor_df$neg_log10_padj <- -log10(pmax(cor_df$padj, .Machine$double.xmin))
  cor_df$significant_BH <- cor_df$padj < 0.05
  cor_df <- cor_df[order(cor_df$rho), ]
  labels <- gsub("_", " ", cor_df$Module, fixed = TRUE)
  cor_df$ModuleLabel <- factor(labels, levels = labels)
  max_abs_rho <- max(abs(cor_df$rho))
  size_label <- "-log10(padj)"

  p_bubble <- ggplot(cor_df, aes(x = 1, y = ModuleLabel)) +
    geom_point(aes(size = neg_log10_padj, fill = rho),
               shape = 21, color = "black", stroke = 0.6, alpha = 0.9) +
    scale_size_continuous(range = c(4, 14), name = size_label) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-max_abs_rho, max_abs_rho),
      name = "Spearman rho"
    ) +
    xlim(0.5, 1.5) +
    labs(title = "AUTS2 Correlation with 15 Mechano-pathways",
         subtitle = "Spearman correlation | BH-adjusted P across 15 pathways",
         x = NULL, y = NULL) +
    guides(fill = guide_colorbar(order = 1, barwidth = 0.6, barheight = 12),
           size = guide_legend(order = 2, override.aes = list(fill = "grey50"))) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey40"),
          axis.text.y = element_text(size = 10, colour = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 1),
          legend.position = "right",
          legend.box = "vertical",
          plot.margin = margin(10, 15, 5, 5))

  if (!dir.exists(output_dir)) dir.create(output_dir)
  if (!dir.exists(output_dir)) stop("Could not create output folder: ", output_dir)
  ggsave(file.path(output_dir, paste0(output_name, ".pdf")), p_bubble, width = 6, height = 7)
  ggsave(file.path(output_dir, paste0(output_name, ".png")), p_bubble,
         width = 6.5, height = 6.5, dpi = 600, bg = "white")
  write.csv(cor_df[, c("Module", "ModuleSize", "rho", "p_val", "padj",
                       "neg_log10_padj", "significant_BH")],
            file.path(output_dir, paste0(output_name, "_Data.csv")), row.names = FALSE)
  message("Pathways with BH-adjusted P < 0.05: ", sum(cor_df$significant_BH), "/15")
  message("Saved in: ", output_dir)
})

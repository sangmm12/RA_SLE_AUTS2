# ========================================================================
# Integrated Enrichment Analysis Script
# Input: Intersection_All.txt (gene symbols)
# Output: Barplot in OE29.Enrich.R style (GO BP/CC/MF + KEGG)
# ========================================================================

# ========================================================================
# Install required packages (uncomment if first run)
# ========================================================================
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("RColorBrewer")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("clusterProfiler")
# BiocManager::install("enrichplot")

# ========================================================================
# Load libraries
# ========================================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# ========================================================================
# Set working directory (change to your own path if needed)
# ========================================================================
setwd("")

# ========================================================================
# Parameter settings
# ========================================================================
pvalueFilter <- 0.05
p.adjustFilter <- 1
showCategory <- 5  # Number of top terms per ontology to display

# ========================================================================
# Load intersection genes
# ========================================================================
message("Loading IntersectionGenes.txt ...")
gene_symbols <- read.table("IntersectionGenes.txt", header = FALSE,
                           stringsAsFactors = FALSE)[, 1]
gene_symbols <- trimws(gene_symbols)
gene_symbols <- unique(gene_symbols[gene_symbols != ""])
message(sprintf("Total unique genes: %d", length(gene_symbols)))
print(gene_symbols)

# ========================================================================
# Convert gene symbols to Entrez IDs
# ========================================================================
message("Converting gene symbols to Entrez IDs ...")
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene <- unique(entrez_ids$ENTREZID)
message(sprintf("Successfully converted: %d genes", length(gene)))

# ========================================================================
# GO Enrichment Analysis
# ========================================================================
message("Running GO enrichment analysis ...")
kk_go <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db,
                  pvalueCutoff = 1, qvalueCutoff = 1,
                  ont = "all", readable = TRUE)
GO <- as.data.frame(kk_go)
GO <- GO[(GO$pvalue < pvalueFilter & GO$p.adjust < p.adjustFilter), ]
message(sprintf("GO significant terms: %d", nrow(GO)))

# Save GO results
write.table(GO, file = "GO_Enrichment_Results.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: GO_Enrichment_Results.txt")

# ========================================================================
# KEGG Enrichment Analysis
# ========================================================================
message("Running KEGG enrichment analysis ...")
options(timeout = 600)
kk_kegg <- enrichKEGG(gene = gene, organism = "hsa",
                      pvalueCutoff = 1, qvalueCutoff = 1)
KEGG <- as.data.frame(kk_kegg)
KEGG <- KEGG[(KEGG$pvalue < pvalueFilter & KEGG$p.adjust < p.adjustFilter), ]

# Convert KEGG gene IDs back to gene symbols
if (nrow(KEGG) > 0) {
  # Map Entrez IDs in geneID column to gene symbols
  id_map <- setNames(entrez_ids$SYMBOL, entrez_ids$ENTREZID)
  KEGG$geneID <- sapply(strsplit(KEGG$geneID, "/"), function(x) {
    symbols <- id_map[x]
    symbols[is.na(symbols)] <- x[is.na(symbols)]
    paste(symbols, collapse = "/")
  })
}
message(sprintf("KEGG significant terms: %d", nrow(KEGG)))

# Save KEGG results
write.table(KEGG, file = "KEGG_Enrichment_Results.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: KEGG_Enrichment_Results.txt")

# ========================================================================
# Combine GO and KEGG for visualization
# ========================================================================
if (nrow(KEGG) > 0) {
  KEGG_plot <- KEGG[, -c(1, 2)]  # Remove ID and Description? No, keep Description
  KEGG_plot <- cbind(ONTOLOGY = "KEGG", KEGG_plot)
  GK <- rbind(GO, KEGG_plot)
} else {
  GK <- GO
  GK$ONTOLOGY <- as.character(GK$ONTOLOGY)
}

GK$ONTOLOGY <- factor(GK$ONTOLOGY, levels = c("BP", "CC", "MF", "KEGG"))

# ========================================================================
# Prepare data for barplot (OE29 style)
# ========================================================================
# Select top N terms per ontology by p-value, remove NAs
data_selected <- GK %>%
  filter(!is.na(ONTOLOGY)) %>%
  arrange(ONTOLOGY, pvalue) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = showCategory) %>%
  ungroup()

# Abbreviate long description labels
data_selected$Description_short <- data_selected$Description
# Truncate descriptions longer than 60 characters for cleaner display
long_idx <- nchar(data_selected$Description_short) > 60
data_selected$Description_short[long_idx] <- paste0(
  substr(data_selected$Description_short[long_idx], 1, 57), "...")

data_selected$log10p <- -log10(data_selected$pvalue)
data_selected <- data_selected %>%
  arrange(ONTOLOGY, desc(log10p)) %>%
  mutate(Description_short = factor(Description_short,
                                    levels = rev(unique(Description_short))))

message(sprintf("Total terms for plotting: %d", nrow(data_selected)))

# ========================================================================
# Define ontology colors (matching OE29.Enrich.R style)
# ========================================================================
ontology_colors <- c("BP" = "#5B8FA8",
                     "CC" = "#82B69A",
                     "MF" = "#B49FCC",
                     "KEGG" = "#E5A06B")

# Ensure all levels present
used_levels <- levels(data_selected$ONTOLOGY)
plot_colors <- ontology_colors[names(ontology_colors) %in% used_levels]

# ========================================================================
# Calculate group positions for side labels
# ========================================================================
group_positions <- data_selected %>%
  group_by(ONTOLOGY) %>%
  summarise(
    y_start = min(as.numeric(Description_short)),
    y_end = max(as.numeric(Description_short)),
    .groups = "drop"
  ) %>%
  mutate(
    y_start = y_start - 0.4,
    y_end = y_end + 0.4
  )

# ========================================================================
# Create barplot (OE29.Enrich.R style)
# ========================================================================
message("Creating barplot ...")

# Calculate x-axis range
max_log10p <- max(data_selected$log10p, na.rm = TRUE)
x_left_limit <- -max_log10p * 0.20

p <- ggplot(data_selected, aes(x = log10p, y = Description_short, fill = ONTOLOGY)) +
  geom_col(width = 0.8) +
  # Colored side panels
  geom_rect(data = group_positions,
            aes(xmin = -max_log10p * 0.15, xmax = -max_log10p * 0.07,
                ymin = y_start, ymax = y_end, fill = ONTOLOGY),
            inherit.aes = FALSE, alpha = 0.8) +
  # Ontology labels on side panels
  geom_text(data = group_positions,
            aes(x = -max_log10p * 0.11, y = (y_start + y_end) / 2,
                label = ONTOLOGY),
            inherit.aes = FALSE, size = 4, fontface = "bold", color = "black") +
  # Description labels inside bars
  geom_text(aes(label = Description_short, x = 0.05),
            hjust = 0, size = 3.5, color = "black", fontface = "plain") +
  scale_fill_manual(values = plot_colors, name = "Ontology") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.05)),
                     limits = c(x_left_limit, NA)) +
  labs(title = "Enrichment Analysis (Intersection Genes)",
       x = "-log10(pvalue)", y = "") +
  theme_void(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt"),
    axis.line.x = element_line(color = "black", size = 0.3),
    axis.ticks.x = element_line(color = "black"),
    axis.text.x = element_text(size = 10)
  )

# Save barplot
n_terms <- nrow(data_selected)
plot_height <- max(6, n_terms * 0.45)

pdf(file = "Enrichment_Barplot.pdf", width = 9, height = plot_height)
print(p)
dev.off()
message(sprintf("Saved: Enrichment_Barplot.pdf (%.1f x %.1f)", 9, plot_height))

# ========================================================================
# Also save as PNG for quick preview
# ========================================================================
tryCatch({
  png(file = "Enrichment_Barplot.png", width = 9, height = plot_height,
      units = "in", res = 300)
  print(p)
  dev.off()
  message("Saved: Enrichment_Barplot.png")
}, error = function(e) {
  message(sprintf("PNG save skipped: %s", e$message))
})

# ========================================================================
# Summary
# ========================================================================
message("\n========== ANALYSIS SUMMARY ==========")
message(sprintf("Input genes: %d", length(gene_symbols)))
message(sprintf("Converted genes: %d", length(gene)))
message(sprintf("GO terms: %d", nrow(GO)))
if (nrow(KEGG) > 0) {
  message(sprintf("KEGG pathways: %d", nrow(KEGG)))
} else {
  message("KEGG pathways: 0 (none significant)")
}
message(sprintf("Terms in barplot: %d", n_terms))
message("\nAll analysis complete!")


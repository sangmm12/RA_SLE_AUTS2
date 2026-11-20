

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(WGCNA)
  library(ggnewscale)
})
if (!requireNamespace("hdWGCNA", quietly=TRUE)) stop("hdWGCNA is required")
# 提供原脚本使用的 hdWGCNA 函数名，避免 library(hdWGCNA) 在旧版本中
# 附带加载 enrichR 并访问远程服务器。函数本身仍来自已安装的 hdWGCNA。
for (nm in getNamespaceExports("hdWGCNA")) {
  assign(nm, getExportedValue("hdWGCNA", nm), envir=environment())
}


# ============================================================================
# 
# ============================================================================
cat_construct_metacells <- function(seurat_obj, k, name, min_cells = 50, max_shared = 10, assay = "RNA") {
  if ("cell_type" %in% colnames(seurat_obj@meta.data) &
      "donor" %in% colnames(seurat_obj@meta.data)) {
    set.seed(717)
    seurat_obj <- hdWGCNA::SetupForWGCNA(
      seurat_obj,
      gene_select = "fraction",
      fraction = 0.05,
      wgcna_name = name
    )
    seurat_obj <- hdWGCNA::MetacellsByGroups(
      seurat_obj = seurat_obj,
      group.by = c("cell_type", "donor"),
      k = k,
      min_cells = min_cells,
      max_shared = max_shared,
      mode = "sum",
      assay = assay,
      ident.group = "cell_type"
    )
    seurat_obj <- hdWGCNA::NormalizeMetacells(seurat_obj)
    return(seurat_obj)
  }
  stop("The column name of cell_type or donor does not exist in the meta.data of your Seurat object, please add it!")
}

ConstructNetwork <- function (seurat_obj, soft_power = NULL, tom_outdir = "TOM", 
                              use_metacells = TRUE, setDatExpr = TRUE, group.by = NULL, 
                              group_name = NULL, consensus = FALSE, multi.group.by = NULL, 
                              multi_groups = NULL, blocks = NULL, maxBlockSize = 30000, 
                              randomSeed = 12345, corType = "pearson", consensusQuantile = 0.3, 
                              networkType = "signed", TOMType = "signed", TOMDenom = "min", 
                              scaleTOMs = TRUE, scaleQuantile = 0.8, sampleForScaling = TRUE, 
                              sampleForScalingFactor = 1000, useDiskCache = TRUE, chunkSize = NULL, 
                              deepSplit = 4, pamStage = FALSE, detectCutHeight = 0.995, 
                              minModuleSize = 50, mergeCutHeight = 0.2, saveConsensusTOMs = TRUE, nThreads = 6,
                              ...) 
{
  if (consensus) {
    if (!("multiExpr" %in% names(GetActiveWGCNA(seurat_obj))) | 
        setDatExpr == TRUE) {
      seurat_obj <- SetMultiExpr(seurat_obj, group_name = group_name, 
                                 group.by = group.by, multi.group.by = multi.group.by, 
                                 multi_groups = multi_groups)
    }
    multiExpr <- GetMultiExpr(seurat_obj)
    checkSets(multiExpr)
  }
  else {
    if (!("datExpr" %in% names(GetActiveWGCNA(seurat_obj))) | 
        setDatExpr == TRUE) {
      seurat_obj <- SetDatExpr(seurat_obj, group_name = group_name, 
                               group.by = group.by, use_metacells = use_metacells, 
                               return_seurat = TRUE)
    }
    datExpr <- GetDatExpr(seurat_obj)
    if (is.null(group_name)) {
      group_name <- "all"
    }
    nSets = 1
    setLabels = gsub(" ", "_", group_name)
    shortLabels = setLabels
    multiExpr <- list()
    multiExpr[[group_name]] <- list(data = datExpr)
    checkSets(multiExpr)
  }
  if (!dir.exists(tom_outdir)) {
    dir.create(tom_outdir)
  }
  if (is.null(soft_power) & !consensus) {
    soft_power <- GetPowerTable(seurat_obj) %>% subset(SFT.R.sq >= 
                                                         0.8) %>% .$Power %>% min
    cat(paste0("Soft power not provided. Automatically using the lowest power that meets 0.8 scale-free topology fit. Using soft_power = ", 
               soft_power, "\n"))
  }
  else if (consensus) {
    power_tables <- GetPowerTable(seurat_obj) %>% dplyr::group_split(group)
    soft_power <- sapply(power_tables, function(power_table) {
      power_table %>% subset(SFT.R.sq >= 0.8) %>% .$Power %>% 
        min
    })
    cat(paste0("Soft power not provided. Automatically using the lowest power that meets 0.8 scale-free topology fit. Using soft_power = c(", 
               paste0(soft_power, collapse = ","), ")\n"))
  }
  net <- WGCNA::blockwiseConsensusModules(multiExpr, power = soft_power, 
                                          blocks = blocks, maxBlockSize = maxBlockSize, randomSeed = randomSeed, 
                                          corType = corType, consensusQuantile = consensusQuantile, 
                                          networkType = networkType, TOMType = TOMType, TOMDenom = TOMDenom, 
                                          scaleTOMs = scaleTOMs, scaleQuantile = scaleQuantile, 
                                          sampleForScaling = sampleForScaling, sampleForScalingFactor = sampleForScalingFactor, 
                                          useDiskCache = useDiskCache, chunkSize = chunkSize, 
                                          deepSplit = deepSplit, pamStage = pamStage, detectCutHeight = detectCutHeight, 
                                          minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight, 
                                          saveConsensusTOMs = saveConsensusTOMs, consensusTOMFilePattern = "ConsensusTOM-block.%b.rda", 
                                          ...)
  file.rename("ConsensusTOM-block.1.rda", paste0(tom_outdir, "/", gsub(" ", 
                                                                       "_", group_name), "_ConsensusTOM-block.1.rda"))
  params <- list(power = soft_power, blocks = blocks, maxBlockSize = maxBlockSize, 
                 randomSeed = randomSeed, corType = corType, consensusQuantile = consensusQuantile, 
                 networkType = networkType, TOMType = TOMType, TOMDenom = TOMDenom, 
                 scaleTOMs = scaleTOMs, scaleQuantile = scaleQuantile, 
                 sampleForScaling = sampleForScaling, sampleForScalingFactor = sampleForScalingFactor, 
                 useDiskCache = useDiskCache, chunkSize = chunkSize, 
                 deepSplit = deepSplit, pamStage = pamStage, detectCutHeight = detectCutHeight, 
                 minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight, 
                 saveConsensusTOMs = saveConsensusTOMs)
  seurat_obj <- SetWGCNAParams(seurat_obj, params)
  net$TOMFiles <- paste0(getwd(), '/', tom_outdir, '/', group_name, '_', net$TOMFiles)
  
  seurat_obj <- SetNetworkData(seurat_obj, net)
  mods <- GetNetworkData(seurat_obj)$colors
  seurat_obj <- SetModules(seurat_obj, modules = data.frame(gene_name = names(mods), 
                                                           module = factor(mods, levels = unique(mods)), color = mods))
  seurat_obj
}


outdir <- file.path(result_dir, "F_network")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
if (run_network_from_raw) {
  seurat_obj_path <- file.path(root_dir, "7pbmc_celltype.rds")
  k <- 20
# ---- 2. 读取数据并提取成纤维细胞 ----
cat(">>> 读取数据...\n")
seurat_obj <- readRDS(seurat_obj_path)

cat(">>> 数据中的细胞类型:\n")
print(table(seurat_obj$celltype.main))

fibro_type <- "Fibroblast"
stopifnot("未找到 Fibroblast 细胞类型" = fibro_type %in% seurat_obj$celltype.main)

# ---- 3. 提取成纤维细胞子集 ----
cat(">>> 提取成纤维细胞子集...\n")
fibro_obj <- subset(seurat_obj, celltype.main == fibro_type)
cat(">>> 成纤维细胞数量:", ncol(fibro_obj), "\n")

# ---- 4. 添加 donor 列 ----
fibro_obj$cell_type <- fibro_obj$celltype.main
fibro_obj$donor <- fibro_obj$sample
cat(">>> donor 分布:\n")
print(table(fibro_obj$donor))

# ---- 5. 构建Metacells ----
cat(">>> 构建 Metacells (k =", k, ")...\n")
fibro_obj <- cat_construct_metacells(
  seurat_obj = fibro_obj,
  k = k,
  name = "Fibroblast"
)

metacell_obj <- GetMetacellObject(fibro_obj)
cat(">>> Metacell 对象维度:", dim(metacell_obj), "\n")
cat(">>> 各 donor 的 metacell 数量:\n")
print(table(metacell_obj$donor))

saveRDS(fibro_obj, file.path(outdir, "fibroblast_metacell.rds"))

# ---- 6. WGCNA: 设置表达数据 ----
cat(">>> 设置表达数据...\n")
fibro_obj <- SetDatExpr(
  fibro_obj,
  group_name = fibro_type,
  group.by = "cell_type",
  assay = "RNA"
)

# ---- 7. 软阈值选择 ----
cat(">>> 测试软阈值...\n")
fibro_obj <- TestSoftPowers(fibro_obj)

plot_list <- PlotSoftPowers(fibro_obj)
p_soft <- wrap_plots(plot_list, ncol = 2)
ggsave(file.path(outdir, "01_soft_power.pdf"),
       plot = p_soft, height = 6, width = 10)

power_table <- GetPowerTable(fibro_obj)
write_csv(power_table, file.path(outdir, "soft_power_table.csv"))
cat(">>> 软阈值表:\n")
print(head(power_table))

# ---- 8. 构建网络 ----
cat(">>> 构建共表达网络...\n")
auto_power <- power_table %>%
  subset(SFT.R.sq >= 0.8) %>%
  .$Power %>%
  min()
cat(">>> 自动选择 soft_power =", auto_power, "\n")

fibro_obj <- ConstructNetwork(
  fibro_obj,
  soft_power = auto_power,
  setDatExpr = FALSE,
  tom_outdir = outdir,
  nThreads = 10
)

pdf(file.path(outdir, "02_dendrogram.pdf"), width = 8, height = 5)
PlotDendrogram(fibro_obj, main = "Fibroblast hdWGCNA Dendrogram")
dev.off()

# TOM矩阵加载（直接load实际文件，自动获取变量名）
tom_file <- file.path(outdir, "all_ConsensusTOM-block.1.rda")
if (file.exists(tom_file)) {
  cat(">>> 找到TOM文件:", tom_file, "\n")
  tom_vars <- load(tom_file)
  cat(">>> TOM文件中的变量:", tom_vars, "\n")
  tom_obj <- get(tom_vars[1])
  TOM <- as.matrix(tom_obj)
} else {
  TOM <- GetTOM(fibro_obj)
}

# ---- 9. 模块特征基因 ----
cat(">>> 计算模块特征基因...\n")
fibro_obj <- ModuleEigengenes(
  fibro_obj,
  group.by.vars = "donor"
)

hMEs <- GetMEs(fibro_obj, harmonized = TRUE)
# 排除grey模块
if ("grey" %in% colnames(hMEs)) hMEs <- hMEs[, colnames(hMEs) != "grey"]
cat(">>> 参与分析的模块:", paste(colnames(hMEs), collapse = ", "), "\n")
MEs  <- GetMEs(fibro_obj, harmonized = FALSE)

# ---- 10. 模块连通性 ----
cat(">>> 计算模块连通性...\n")
fibro_obj <- ModuleConnectivity(
  fibro_obj,
  group.by = "cell_type",
  group_name = fibro_type
)

p_kme <- PlotKMEs(fibro_obj, ncol = 5)
ggsave(file.path(outdir, "03_KME_plots.pdf"),
       plot = p_kme, height = 6, width = 12)

# ---- 11. 导出模块基因表 ----
modules <- GetModules(fibro_obj)
write_csv(modules, file.path(outdir, "modules.csv"))
cat(">>> 模块组成:\n")
print(table(modules$module))


  saveRDS(fibro_obj, file.path(outdir, "fibroblast_wgcna.rds"))
  input_wgcna_dir <- outdir
} else {
  fibro_obj <- readRDS(file.path(input_wgcna_dir, "fibroblast_wgcna.rds"))
  pdf(file.path(outdir, "02_dendrogram.pdf"), width=8, height=5)
  PlotDendrogram(fibro_obj, main="Fibroblast hdWGCNA Dendrogram")
  dev.off()
}
wgcna_obj_path <- file.path(input_wgcna_dir, "fibroblast_wgcna.rds")

# ============================================================================
#
# ============================================================================
(function() {
# ============================================================================
===========================================================================

options(stringsAsFactors = FALSE)

# ---- 0. User parameters -----------------------------------------------------
target_gene <- "AUTS2"

# The script is saved in the project root. If you run it from another folder,
# leave project_dir as NULL so the script path is detected automatically.
project_dir <- root_dir
wgcna_dir_name <- "Fibroblast_WGCNA"

assay_name <- "RNA"
assay_layer <- "data"
exclude_grey_module <- TRUE

# "both": p < threshold and abs(r) >= threshold
# "positive": p < threshold and r >= threshold
# "negative": p < threshold and r <= -threshold
coexpression_mode <- "both"
cor_method <- "pearson"
r_threshold <- 0.30
p_threshold <- 0.05
use_fdr <- FALSE

# Set TRUE to force recalculation from the Seurat metacell object.
# Set FALSE to reuse AUTS2_all_gene_correlation.csv when available.
refresh_correlation <- FALSE

figure_basename <- "14_AUTS2_module_coexpression_bubble"
figure_width_mm <- 56
figure_height_mm <- 68
figure_dpi <- 600
x_axis_label <- "Smooth muscle cell"
size_legend_title <- "Percent expressed"
color_legend_title <- "Average expression"

# WGCNA-like order, shown from top to bottom when present.
module_order_top_to_bottom <- c(
  "tan", "pink", "green", "greenyellow", "black", "brown", "cyan",
  "magenta", "red", "yellow", "purple", "salmon", "turquoise", "blue"
)

# ---- 1. Helpers -------------------------------------------------------------
message_step <- function(...) {
  message(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = ""))
}

stop_missing_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing R package(s): ", paste(missing, collapse = ", "), "\n",
      "Install with:\ninstall.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "),
      "))",
      call. = FALSE
    )
  }
}

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  }

  frame_files <- vapply(sys.frames(), function(x) {
    if (!is.null(x$ofile)) x$ofile else NA_character_
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    return(dirname(normalizePath(frame_files[1], winslash = "/")))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (!is.null(active_path) && nzchar(active_path)) {
      return(dirname(normalizePath(active_path, winslash = "/")))
    }
  }

  normalizePath(getwd(), winslash = "/")
}

theme_pub_dot <- function(base_size = 7, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = base_size, colour = "black"),
      axis.text.y = ggplot2::element_text(size = base_size, colour = "black"),
      legend.title = ggplot2::element_text(size = base_size - 0.3),
      legend.text = ggplot2::element_text(size = base_size - 0.7),
      legend.key.height = grid::unit(3.8, "mm"),
      legend.key.width = grid::unit(3.8, "mm"),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.35),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 3, 3, 3, unit = "mm")
    )
}

save_figure <- function(plot, basename, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  ggplot2::ggsave(
    filename = paste0(basename, ".pdf"),
    plot = plot,
    width = width_in,
    height = height_in,
    device = grDevices::cairo_pdf,
    family = "Arial",
    units = "in"
  )

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(
      filename = paste0(basename, ".svg"),
      plot = plot,
      width = width_in,
      height = height_in,
      device = svglite::svglite,
      units = "in"
    )
  } else {
    warning("Package 'svglite' is not installed; SVG export skipped.")
  }

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(
      filename = paste0(basename, ".tiff"),
      plot = plot,
      width = width_in,
      height = height_in,
      device = ragg::agg_tiff,
      units = "in",
      dpi = dpi,
      compression = "lzw"
    )
    ggplot2::ggsave(
      filename = paste0(basename, ".png"),
      plot = plot,
      width = width_in,
      height = height_in,
      device = ragg::agg_png,
      units = "in",
      dpi = dpi
    )
  } else {
    warning("Package 'ragg' is not installed; TIFF/PNG export skipped.")
  }
}

read_expression_matrix <- function(metacell_obj, assay = "RNA", layer = "data") {
  tryCatch(
    Seurat::GetAssayData(metacell_obj, assay = assay, layer = layer),
    error = function(e_layer) {
      Seurat::GetAssayData(metacell_obj, assay = assay, slot = layer)
    }
  )
}

match_gene_name <- function(gene, row_names) {
  hit <- row_names[row_names == gene]
  if (length(hit) > 0) return(hit[1])

  hit <- row_names[toupper(row_names) == toupper(gene)]
  if (length(hit) > 0) return(hit[1])

  stop("Target gene not found in expression matrix: ", gene, call. = FALSE)
}

compute_target_correlations <- function(expr_matrix, modules, target_gene,
                                        method = "pearson") {
  target_in_matrix <- match_gene_name(target_gene, rownames(expr_matrix))

  if (exclude_grey_module) {
    modules <- modules |> dplyr::filter(.data$module != "grey")
  }
  module_genes <- modules |>
    dplyr::pull("gene_name") |>
    unique()
  module_genes <- setdiff(module_genes, target_in_matrix)
  module_genes <- intersect(module_genes, rownames(expr_matrix))

  if (length(module_genes) == 0) {
    stop("No module genes were found in the expression matrix.", call. = FALSE)
  }

  target_expr <- as.numeric(expr_matrix[target_in_matrix, ])
  gene_expr <- expr_matrix[module_genes, , drop = FALSE]

  message_step("Computing ", length(module_genes), " gene correlations with ", target_in_matrix)

  cor_list <- lapply(module_genes, function(gene) {
    y <- as.numeric(gene_expr[gene, ])
    ok <- stats::complete.cases(target_expr, y)
    n_ok <- sum(ok)
    if (n_ok < 5 || stats::sd(y[ok]) == 0 || stats::sd(target_expr[ok]) == 0) {
      return(tibble::tibble(gene = gene, estimate = NA_real_, p_value = NA_real_, n = n_ok))
    }

    test <- tryCatch(
      stats::cor.test(target_expr[ok], y[ok], method = method),
      error = function(e) NULL
    )
    if (is.null(test)) {
      return(tibble::tibble(gene = gene, estimate = NA_real_, p_value = NA_real_, n = n_ok))
    }

    tibble::tibble(
      gene = gene,
      estimate = unname(test$estimate),
      p_value = test$p.value,
      n = n_ok
    )
  })

  dplyr::bind_rows(cor_list) |>
    dplyr::mutate(fdr = stats::p.adjust(.data$p_value, method = "BH")) |>
    dplyr::arrange(.data$p_value)
}

flag_coexpressed <- function(df, mode = "both", r_cutoff = 0.30,
                             p_cutoff = 0.05, adjusted = FALSE) {
  p_col <- if (adjusted) "fdr" else "p_value"
  if (!p_col %in% colnames(df)) {
    stop("Column not found for significance testing: ", p_col, call. = FALSE)
  }

  df |>
    dplyr::mutate(
      is_significant = .data[[p_col]] < p_cutoff,
      is_coexpressed = dplyr::case_when(
        mode == "both" ~ .data$is_significant & abs(.data$estimate) >= r_cutoff,
        mode == "positive" ~ .data$is_significant & .data$estimate >= r_cutoff,
        mode == "negative" ~ .data$is_significant & .data$estimate <= -r_cutoff,
        TRUE ~ NA
      ),
      direction = dplyr::case_when(
        .data$is_coexpressed & .data$estimate > 0 ~ "positive",
        .data$is_coexpressed & .data$estimate < 0 ~ "negative",
        TRUE ~ "not_selected"
      )
    )
}

# ---- 2. Load packages and paths --------------------------------------------
stop_missing_packages(c("ggplot2", "dplyr", "readr", "tibble", "scales"))

library(ggplot2)
library(dplyr)
library(readr)
library(tibble)

if (is.null(project_dir)) {
  project_dir <- get_script_dir()
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
wgcna_dir <- input_wgcna_dir
g_output_dir <- file.path(result_dir, "G_original_bubble")
dir.create(g_output_dir, recursive=TRUE, showWarnings=FALSE)

if (!dir.exists(wgcna_dir)) {
  stop("WGCNA output directory not found: ", wgcna_dir, call. = FALSE)
}

modules_path <- file.path(wgcna_dir, "modules.csv")
corr_path <- file.path(wgcna_dir, paste0(target_gene, "_all_gene_correlation.csv"))
summary_path <- file.path(g_output_dir, paste0(target_gene, "_module_coexpression_bubble_source.csv"))
gene_source_path <- file.path(g_output_dir, paste0(target_gene, "_gene_correlation_with_module.csv"))

if (!file.exists(modules_path)) {
  stop("Module table not found: ", modules_path, call. = FALSE)
}

message_step("Project directory: ", project_dir)
message_step("WGCNA directory: ", wgcna_dir)

modules <- readr::read_csv(modules_path, show_col_types = FALSE) |>
  dplyr::mutate(module = as.character(.data$module))

required_module_cols <- c("gene_name", "module")
if (!all(required_module_cols %in% colnames(modules))) {
  stop("modules.csv must contain columns: gene_name, module", call. = FALSE)
}

if (exclude_grey_module) {
  modules_for_plot <- modules |> dplyr::filter(.data$module != "grey")
} else {
  modules_for_plot <- modules
}

# ---- 3. Get or compute AUTS2-gene correlations -----------------------------
if (!refresh_correlation && file.exists(corr_path)) {
  message_step("Reading existing correlation table: ", basename(corr_path))
  gene_cor <- readr::read_csv(corr_path, show_col_types = FALSE)

  if (!"fdr" %in% colnames(gene_cor)) {
    gene_cor <- gene_cor |>
      dplyr::mutate(fdr = stats::p.adjust(.data$p_value, method = "BH"))
  }
} else {
  stop_missing_packages("Seurat")

  metacell_path <- file.path(wgcna_dir, "fibroblast_metacell.rds")
  wgcna_obj_path <- file.path(wgcna_dir, "fibroblast_wgcna.rds")

  if (file.exists(metacell_path)) {
    message_step("Reading metacell object: ", basename(metacell_path))
    metacell_obj <- readRDS(metacell_path)
  } else if (file.exists(wgcna_obj_path)) {
    message_step("Reading WGCNA object and extracting metacell object.")
    fibro_obj <- readRDS(wgcna_obj_path)
    if (!exists("GetMetacellObject", mode = "function")) {
      custom_function_path <- file.path(project_dir, "custom_function.R")
      if (file.exists(custom_function_path)) {
        source(custom_function_path)
      }
    }
    if (!exists("GetMetacellObject", mode = "function")) {
      stop("GetMetacellObject() is not available and fibroblast_metacell.rds was not found.", call. = FALSE)
    }
    metacell_obj <- GetMetacellObject(fibro_obj)
  } else {
    stop("Neither fibroblast_metacell.rds nor fibroblast_wgcna.rds was found.", call. = FALSE)
  }

  expr_matrix <- read_expression_matrix(metacell_obj, assay = assay_name, layer = assay_layer)
  gene_cor <- compute_target_correlations(expr_matrix, modules_for_plot, target_gene, method = cor_method)
  readr::write_csv(gene_cor, file.path(g_output_dir, basename(corr_path)))
  message_step("Saved correlation table: ", corr_path)
}

required_cor_cols <- c("gene", "estimate", "p_value")
if (!all(required_cor_cols %in% colnames(gene_cor))) {
  stop("Correlation table must contain columns: gene, estimate, p_value", call. = FALSE)
}

gene_cor <- gene_cor |>
  dplyr::left_join(
    modules |> dplyr::select(dplyr::all_of(c("gene_name", "module"))),
    by = c("gene" = "gene_name")
  ) |>
  dplyr::filter(!is.na(.data$module))

if (exclude_grey_module) {
  gene_cor <- gene_cor |> dplyr::filter(.data$module != "grey")
}

gene_cor <- flag_coexpressed(
  gene_cor,
  mode = coexpression_mode,
  r_cutoff = r_threshold,
  p_cutoff = p_threshold,
  adjusted = use_fdr
)

readr::write_csv(gene_cor, gene_source_path)
message_step("Saved gene-level source table: ", gene_source_path)

# ---- 4. Summarise modules ---------------------------------------------------
module_summary <- gene_cor |>
  dplyr::group_by(.data$module) |>
  dplyr::summarise(
    n_genes = dplyr::n(),
    n_coexpressed = sum(.data$is_coexpressed, na.rm = TRUE),
    n_positive = sum(.data$direction == "positive", na.rm = TRUE),
    n_negative = sum(.data$direction == "negative", na.rm = TRUE),
    mean_r_all = mean(.data$estimate, na.rm = TRUE),
    mean_abs_r_all = mean(abs(.data$estimate), na.rm = TRUE),
    mean_r_coexpressed = if (
      sum(.data$is_coexpressed, na.rm = TRUE) > 0
    ) {
      mean(estimate[is_coexpressed %in% TRUE], na.rm = TRUE)
    } else {
      NA_real_
    },
    mean_abs_r_coexpressed = if (
      sum(.data$is_coexpressed, na.rm = TRUE) > 0
    ) {
      mean(abs(estimate[is_coexpressed %in% TRUE]), na.rm = TRUE)
    } else {
      0
    },
    .groups = "drop"
  ) |>
  dplyr::mutate(
    percent_coexpressed = 100 * .data$n_coexpressed / .data$n_genes
  ) |>
  dplyr::arrange(dplyr::desc(.data$percent_coexpressed), .data$module)

readr::write_csv(module_summary, summary_path)
message_step("Saved module-level source table: ", summary_path)

# ---- 5. Build the bubble plot ----------------------------------------------
present_order <- module_order_top_to_bottom[module_order_top_to_bottom %in% module_summary$module]
extra_order <- setdiff(module_summary$module, present_order)
extra_order <- module_summary |>
  dplyr::filter(.data$module %in% extra_order) |>
  dplyr::arrange(dplyr::desc(.data$percent_coexpressed), .data$module) |>
  dplyr::pull("module")

plot_order_top_to_bottom <- c(present_order, extra_order)

plot_df <- module_summary |>
  dplyr::mutate(
    module = factor(.data$module, levels = rev(plot_order_top_to_bottom)),
    target = x_axis_label,
    color_value = .data$mean_abs_r_all
  )

max_percent <- max(plot_df$percent_coexpressed, na.rm = TRUE)
if (!is.finite(max_percent) || max_percent <= 0) max_percent <- 1
size_breaks <- pretty(c(0, max_percent), n = 4)
size_breaks <- size_breaks[size_breaks > 0]
if (length(size_breaks) == 0) size_breaks <- max_percent
max_color <- max(plot_df$color_value, na.rm = TRUE)
if (!is.finite(max_color) || max_color <= 0) max_color <- 1

plot_title <- NULL
size_title <- size_legend_title
color_title <- color_legend_title

p_bubble <- ggplot(
  plot_df,
  aes(x = .data$target, y = .data$module)
) +
  geom_point(
    aes(size = .data$percent_coexpressed, fill = .data$color_value),
    shape = 21,
    colour = "white",
    stroke = 0.18,
    alpha = 0.98
  ) +
  scale_size_area(
    name = size_title,
    max_size = 7,
    breaks = size_breaks,
    labels = function(x) sprintf("%.1f", x)
  ) +
  scale_fill_gradient(
    name = color_title,
    low = "#F2F2F2",
    high = "#E31A1C",
    limits = c(0, max_color),
    labels = function(x) sprintf("%.2f", x)
  ) +
  guides(
    size = guide_legend(
      order = 1,
      override.aes = list(fill = "black", colour = "black", alpha = 1)
    ),
    fill = guide_colorbar(
      order = 2,
      barwidth = grid::unit(3.2, "mm"),
      barheight = grid::unit(13, "mm")
    )
  ) +
  labs(title = plot_title) +
  coord_cartesian(clip = "off") +
  theme_pub_dot(base_size = 7) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = grid::unit(1.5, "mm"),
    plot.title = element_blank()
  )

# ---- 6. Export --------------------------------------------------------------
figure_base_path <- file.path(g_output_dir, figure_basename)
save_figure(
  plot = p_bubble,
  basename = figure_base_path,
  width_mm = figure_width_mm,
  height_mm = figure_height_mm,
  dpi = figure_dpi
)

message_step("Done.")
message_step("Figure basename: ", figure_base_path)
message_step("Bubble size: percent of module genes passing the co-expression cutoff.")
message_step("Bubble color: average absolute Pearson r among co-expressed genes.")
message_step("Cutoff: ", ifelse(use_fdr, "FDR", "p"), " < ", p_threshold,
             ", mode = ", coexpression_mode, ", |r|/r cutoff = ", r_threshold)


})()

# ============================================================================
# 
# ============================================================================

outdir <- file.path(result_dir, "H_module_featureplot")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

plot_list <- ModuleFeaturePlot(fibro_obj)
n_mods <- length(plot_list)
ncol_mod <- min(6, ceiling(sqrt(n_mods)))
nrow_mod <- ceiling(n_mods / ncol_mod)
p_modules <- wrap_plots(plot_list, ncol = ncol_mod)
ggsave(file.path(outdir, "09_module_featureplot.pdf"),
       plot = p_modules,
       height = nrow_mod * 3,
       width = ncol_mod * 3.5,
       limitsize = FALSE)
cat(">>> ModuleFeaturePlot 已保存\n")


# 从原代码生成的列表中单独导出参考图使用的 turquoise 面板，样式不变。
if (is.list(plot_list) && "turquoise" %in% names(plot_list)) {
  ggsave(file.path(outdir, "09b_turquoise_featureplot.pdf"),
         plot_list[["turquoise"]], width=3.5, height=3)
}

# ============================================================================
# 原代码来源：6c-HubGene_black_turquoise.R
# ============================================================================

outdir <- file.path(result_dir, "I_hub_networks")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
wgcna_path <- wgcna_obj_path
net_outdir <- file.path(outdir, "ModuleNetworks")
if (!dir.exists(net_outdir)) dir.create(net_outdir)
# ---- 1. 加载并修正TOM路径 ----
cat(">>> 加载WGCNA结果...\n")
fibro_obj <- readRDS(wgcna_path)

tom_correct <- file.path(input_wgcna_dir, "all_ConsensusTOM-block.1.rda")
if (file.exists(tom_correct)) {
  net_data <- GetNetworkData(fibro_obj)
  net_data$TOMFiles <- tom_correct
  fibro_obj <- SetNetworkData(fibro_obj, net_data)
  cat(">>> TOM路径已修正\n")
}

modules <- GetModules(fibro_obj)



# ---- 3. Turquoise模块环状Hub网络 ----
cat("\n>>> 绘制 turquoise 模块环状网络图...\n")
ModuleNetworkPlot(
  fibro_obj,
  outdir = net_outdir,
  mods = "turquoise",
  n_inner = 15,
  n_outer = 30,
  n_conns = Inf,
  plot_size = c(10, 10),
  vertex.label.cex = 0.8
)

# ---- 4. 导出与AUTS2相关的基因 ----
cat("\n>>> 输出模块中与AUTS2最相关的基因...\n")
gene_cor_path <- file.path(input_wgcna_dir, "AUTS2_all_gene_correlation.csv")
if (file.exists(gene_cor_path)) {
  gene_cor <- read_csv(gene_cor_path)
  for (mod_name in c(, "turquoise")) {
    mod_genes <- modules %>% filter(module == mod_name) %>% pull(gene_name)
    top_auts2 <- gene_cor %>%
      filter(gene %in% mod_genes) %>%
      arrange(desc(abs(estimate))) %>%
      select(gene, estimate, p_value) %>%
      head(30)
    write_csv(top_auts2, file.path(outdir, str_c(mod_name, "_top30_AUTS2_cor.csv")))
    cat(">>>", mod_name, "模块中与AUTS2最相关的基因:\n")
    print(top_auts2, n = 30)
    cat("\n")
  }
}

cat("\n>>> ====== 完成! ======\n")
cat(">>> 输出保存在:", net_outdir, "\n")

# ============================================================================
# 原代码来源：Fibroblast_WGCNA/RA_WGCNA_MMP9_top4_padj_results/RA_WGCNA_MMP9_top4_BH.R
# ============================================================================
# Run with Source in RStudio or with Rscript.
local({
  wgcna_dir <- input_wgcna_dir
  gene_dir <- "D:/mechanical-related/111"
  output_dir <- file.path(result_dir, "J_original_top4_BH")
  # Select the first four rows of the MMP9 bubble plot, ordered by signed Spearman rho.
  bubble_file <- file.path(gene_dir, "RA_MMP9_Mechano_raw_p_results",
                           "RA_MMP9_Mechano_raw_p_Bubble_Data.csv")
  bubble <- read.csv(bubble_file)
  stopifnot(nrow(bubble) == 15L, !anyNA(bubble$Pathway), !anyDuplicated(bubble$Pathway),
            all(is.finite(bubble$rho)), all(is.finite(bubble$p_val)))
  top <- head(bubble[order(-bubble$rho, bubble$Pathway), ], 4L)
  pathway_labels <- top$Pathway

  module_order <- c("magenta", "turquoise", "yellow", "black", "pink", "purple",
                    "red", "greenyellow", "green", "blue", "brown")
  module_colors <- setNames(
    c("#A3519C", "#55BFB8", "#E8E849", "#111111", "#F2AEC3", "#6B4A9C",
      "#ED2627", "#A4CC4A", "#70B943", "#245293", "#922B2D"), module_order
  )

  suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
  })
  object_file <- file.path(wgcna_dir, "fibroblast_wgcna.rds")
  gene_files <- file.path(gene_dir, paste0(pathway_labels, ".txt"))
  inputs <- c(object_file, gene_files)
  missing <- inputs[!file.exists(inputs)]
  if (length(missing)) stop("Missing input files: ", paste(missing, collapse = ", "))

  message("Reading stored WGCNA module eigengenes and RNA expression...")
  fibroblasts <- readRDS(object_file)
  eigengenes <- as.data.frame(hdWGCNA::GetMEs(fibroblasts, harmonized = FALSE))
  available_modules <- setdiff(colnames(eigengenes), "grey")
  if (!setequal(available_modules, module_order)) {
    stop("The saved object does not have the expected 11 non-grey modules.")
  }
  expression <- Seurat::GetAssayData(fibroblasts, assay = "RNA", layer = "data")
  cells <- rownames(eigengenes)
  if (is.null(cells) || anyDuplicated(cells) || !all(cells %in% colnames(expression))) {
    stop("The module eigengenes cannot be matched to RNA expression by cell ID.")
  }
  eigengenes <- eigengenes[, module_order, drop = FALSE]
  expression <- expression[, cells, drop = FALSE]

  # Keep original WGCNA scoring: mean normalized RNA expression of the pathway genes.
  # This is a module–pathway analysis, so full pathway memberships are retained.
  matched_genes <- lapply(gene_files, function(file) {
    genes <- trimws(readLines(file, warn = FALSE, encoding = "UTF-8"))
    intersect(genes[nzchar(genes) & !startsWith(genes, "#")], rownames(expression))
  })
  if (any(lengths(matched_genes) < 3L)) stop("Each pathway must contain at least three measured genes.")
  pathway_scores <- lapply(matched_genes, function(genes) {
    Matrix::colMeans(expression[genes, , drop = FALSE])
  })

  results <- expand.grid(module = module_order, pathway = pathway_labels,
                         stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  estimates <- lapply(seq_len(nrow(results)), function(i) {
    index <- match(results$pathway[i], pathway_labels)
    x <- pathway_scores[[index]]
    y <- eigengenes[[results$module[i]]]
    valid <- is.finite(x) & is.finite(y)
    if (sum(valid) < 5L || sd(x[valid]) == 0 || sd(y[valid]) == 0) {
      stop("Pearson correlation cannot be estimated for ", results$module[i],
           " and ", results$pathway[i])
    }
    test <- cor.test(x[valid], y[valid], method = "pearson")
    log_p <- log(2) + pt(-abs(unname(test$statistic)), df = unname(test$parameter),
                         log.p = TRUE)
    c(n_cells = sum(valid), r = unname(test$estimate), p_value = test$p.value,
      log_p = log_p)
  })
  results <- cbind(results, as.data.frame(do.call(rbind, estimates)))
  results$scored_genes <- lengths(matched_genes)[match(results$pathway, pathway_labels)]
  rank_index <- match(results$pathway, pathway_labels)
  results$MMP9_bubble_rank <- rank_index
  results$MMP9_bubble_spearman_rho <- top$rho[rank_index]
  stopifnot(nrow(results) == 44L, all(is.finite(results$p_value)))
  results$padj <- p.adjust(results$p_value, method = "BH")
  # Retain precise magnitudes for extremely small P values that underflow to zero.
  ordered <- order(results$log_p)
  log_bh <- pmin(0, rev(cummin(rev(results$log_p[ordered] +
                                  log(nrow(results) / seq_len(nrow(results)))))))
  results$neg_log10_padj <- -log_bh[order(ordered)] / log(10)
  results$log_p <- NULL
  results$stars <- ifelse(results$padj < 0.001, "***",
                          ifelse(results$padj < 0.01, "**", ifelse(results$padj < 0.05, "*", "")))

  # Plot every tested pair. Use a symmetric color scale covering the observed r values.
  strip <- ".modules"
  row_order <- c(rev(pathway_labels), strip)
  heat <- results
  heat$module <- factor(heat$module, levels = module_order)
  heat$pathway <- factor(heat$pathway, levels = row_order)
  strip_data <- data.frame(
    module = factor(module_order, levels = module_order),
    pathway = factor(strip, levels = row_order), color = unname(module_colors)
  )
  r_limit <- max(0.1, ceiling(max(abs(results$r)) * 10) / 10)
  plot <- ggplot(heat, aes(x = module, y = pathway)) +
    geom_tile(aes(fill = r), colour = "white", linewidth = 0.3) +
    geom_text(aes(label = stars, colour = abs(r) > r_limit * 0.43),
              fontface = "bold", size = 3.2, show.legend = FALSE) +
    scale_colour_manual(values = c("FALSE" = "grey20", "TRUE" = "white")) +
    scale_fill_gradient2(low = "#2C7BB6", mid = "#F7F7F7", high = "#F46D43",
                         limits = c(-r_limit, r_limit), breaks = c(-r_limit, 0, r_limit),
                         name = "Pearson r") +
    ggnewscale::new_scale_fill() +
    geom_tile(data = strip_data, aes(fill = color), colour = "white", linewidth = 0.3) +
    scale_fill_identity() +
    scale_x_discrete(labels = tools::toTitleCase, expand = c(0, 0)) +
    scale_y_discrete(position = "right", labels = function(x) ifelse(x == strip, "", x),
                     expand = c(0, 0)) +
    labs(x = NULL, y = NULL,
         caption = "Top four pathways by MMP9 Spearman rho. Cell-level Pearson tests; BH across 44 pairs.\n* padj < 0.05   ** padj < 0.01   *** padj < 0.001") +
    theme_classic(base_size = 11) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = "grey20"),
          axis.text.y = element_text(colour = "grey20"),
          legend.key.height = grid::unit(0.5, "cm"),
          plot.caption = element_text(hjust = 0, size = 8),
          plot.margin = margin(8, 12, 8, 12))

  if (!dir.exists(output_dir)) dir.create(output_dir)
  if (!dir.exists(output_dir)) stop("Could not create output folder: ", output_dir)
  prefix <- file.path(output_dir, "RA_WGCNA_MMP9_top4_padj")
  write.csv(results, paste0(prefix, "_statistics.csv"), row.names = FALSE)
  ggsave(paste0(prefix, ".pdf"), plot, width = 8.4, height = 3.1)
  ggsave(paste0(prefix, ".png"), plot, width = 8.4, height = 3.1, dpi = 600, bg = "white")
  message("Cells: ", length(cells), "; tested pairs: ", nrow(results),
          "; raw P < 0.05: ", sum(results$p_value < 0.05),
          "; BH padj < 0.05: ", sum(results$padj < 0.05))
  message("Saved in: ", output_dir)
})

# ============================================================================
# 原代码来源：7-module_pathway_correlation_heatmap.R（可选旧版）
# ============================================================================
if (run_legacy_heatmap) {
(function() {

outdir <- file.path(result_dir, "J_legacy_7_script")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
wgcna_path <- wgcna_obj_path
# ---- 1. 定义4个机械力相关通路基因 ----
pathway_list <- list(
  "ECM_assembly" = scan("D:/mechanical-related/111/Extracellular matrix assembly.txt",
                        what = character(), quiet = TRUE),
  "ECM_organization" = scan("D:/mechanical-related/111/Extracellular matrix organization.txt",
                            what = character(), quiet = TRUE),
  "Actin_cytoskeleton" = scan("D:/mechanical-related/111/Actin cytoskeleton.txt",
                              what = character(), quiet = TRUE),
  "Rho_signaling" = scan("D:/mechanical-related/111/Rho protein signal transduction.txt",
                         what = character(), quiet = TRUE)
)

# ---- 2. 加载数据 ----
cat(">>> 加载WGCNA结果...\n")
fibro_obj <- readRDS(wgcna_path)
metacell_obj <- GetMetacellObject(fibro_obj)

hMEs <- GetMEs(fibro_obj, harmonized = FALSE)
if ("grey" %in% colnames(hMEs)) hMEs <- hMEs[, colnames(hMEs) != "grey"]
cat(">>> 参与分析的模块:", paste(colnames(hMEs), collapse = ", "), "\n")
cat(">>> MEs行数:", nrow(hMEs), " 列数:", ncol(hMEs), "\n")
cat(">>> MEs行名前5:", paste(rownames(hMEs)[1:5], collapse = ", "), "\n")

# ---- 3. 过滤基因并计算通路打分（单细胞水平，与hMEs匹配） ----
expr_matrix_sc <- GetAssayData(fibro_obj, layer = "data", assay = "RNA")
cat(">>> 单细胞表达矩阵维度:", dim(expr_matrix_sc), "\n")

for (nm in names(pathway_list)) {
  genes_in <- pathway_list[[nm]][pathway_list[[nm]] %in% rownames(expr_matrix_sc)]
  cat(">>>", nm, ":", length(genes_in), "/", length(pathway_list[[nm]]), "个基因在数据中\n")
  if (length(genes_in) < 3) next
  pathway_score <- colMeans(expr_matrix_sc[genes_in, , drop = FALSE])
  fibro_obj[[nm]] <- pathway_score
}

avail_traits <- names(pathway_list)[sapply(pathway_list, function(g) sum(g %in% rownames(expr_matrix_sc)) >= 3)]
cat(">>> 可用通路:", paste(avail_traits, collapse = ", "), "\n")

# ---- 4. 合并通路打分与MEs（均在同一水平：单细胞） ----
score_df <- fibro_obj@meta.data[, avail_traits, drop = FALSE]
cat(">>> score_df（单细胞水平）行数:", nrow(score_df), " 前3行名:", paste(rownames(score_df)[1:3], collapse = ", "), "\n")

me_df <- as.data.frame(hMEs)
cat(">>> me_df 行数:", nrow(me_df), " 前3行名:", paste(rownames(me_df)[1:3], collapse = ", "), "\n")

cor_df <- merge(score_df, me_df, by = "row.names") %>%
  column_to_rownames("Row.names")
cat(">>> cor_df 行数:", nrow(cor_df), " 列数:", ncol(cor_df), "\n")
cat(">>> cor_df NA数:", sum(is.na(cor_df)), "\n")

# ---- 5. 计算模块-通路相关性 ----
cat(">>> 计算模块-通路相关性...\n")
module_names <- colnames(hMEs)

cor_list <- list()
for (mod in module_names) {
  for (trait in avail_traits) {
    x <- cor_df[[trait]]
    y <- cor_df[[mod]]
    valid <- complete.cases(x, y)
    if (sum(valid) < 5) {
      cat("  !!", mod, "vs", trait, ": 有效观测", sum(valid), "跳过\n")
      next
    }
    test <- cor.test(x[valid], y[valid], method = "pearson")
    cor_list[[length(cor_list) + 1]] <- tibble(
      module = mod, trait = trait,
      estimate = as.numeric(test$estimate), p_value = test$p.value
    )
  }
}

if (length(cor_list) == 0) {
  stop("没有成功计算任何相关性！请检查数据")
}

cor_res <- bind_rows(cor_list) %>%
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
  mutate(stars = case_when(
    fdr < 0.001 ~ "***", fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*", TRUE ~ ""
  ))

write_csv(cor_res, file.path(outdir, "module_pathway_correlation.csv"))
cat(">>> 相关性结果:\n")
print(cor_res %>% arrange(module, trait))

names(cor_res)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(tibble)
library(ggnewscale)

cat(">>> Advanced compact heatmap with significance stars...\n")

# ---------------------------
# 1. data
# ---------------------------
module_levels <- unique(cor_res$module)
trait_levels  <- unique(cor_res$trait)

bar_row <- ".module_color_bar"

module_levels <- unique(cor_res$module)
trait_levels  <- unique(cor_res$trait)

bar_row <- ".module_color_bar"

plot_df <- cor_res %>%
  select(module, trait, estimate, p_value) %>%
  mutate(
    sig_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    ),
    module = factor(module, levels = module_levels),
    trait  = factor(trait, levels = c(bar_row, rev(trait_levels)))
  )

module_color_df <- tibble(
  module = factor(module_levels, levels = module_levels),
  trait = factor(bar_row, levels = c(bar_row, rev(trait_levels))),
  module_color = str_remove(as.character(module_levels), "^ME") %>% str_to_lower(),
  module_label = str_remove(as.character(module_levels), "^ME") %>% str_to_title()
)

corr_lim <- 0.2

# ---------------------------
# 2. plot
# ---------------------------
p_heatmap <- ggplot() +
  
  geom_tile(
    data = plot_df,
    aes(x = module, y = trait, fill = estimate),
    color = "white",
    linewidth = 0.25
  ) +
  
  geom_text(
    data = plot_df,
    aes(x = module, y = trait, label = sig_label),
    color = "white",
    size = 2.6,
    fontface = "bold"
  ) +
  
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "#F7F7F7",
    high = "#F46D43",
    midpoint = 0,
    limits = c(-corr_lim, corr_lim),
    oob = scales::squish,
    name = NULL,
    breaks = c(-corr_lim, 0, corr_lim),
    labels = c(paste0("-", corr_lim), "0", as.character(corr_lim))
  ) +
  
  ggnewscale::new_scale_fill() +
  
  geom_tile(
    data = module_color_df,
    aes(x = module, y = trait, fill = module_color),
    color = "white",
    linewidth = 0.25
  ) +
  
  scale_fill_identity() +
  
  scale_x_discrete(
    labels = module_color_df$module_label,
    expand = c(0, 0)
  ) +
  
  scale_y_discrete(
    position = "right",
    labels = function(x) ifelse(x == bar_row, "", x),
    expand = c(0, 0)
  ) +
  
  coord_fixed(ratio = 0.75, clip = "off") +
  
  labs(
    title = "",
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(base_size = 10) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    
    axis.text.x = element_text(
      angle = 50,
      hjust = 1,
      vjust = 1,
      size = 8,
      color = "grey25"
    ),
    
    axis.text.y = element_text(
      size = 9,
      color = "grey20",
      margin = margin(l = 4)
    ),
    
    legend.position = "right",
    legend.key.height = unit(0.45, "cm"),
    legend.key.width = unit(0.18, "cm"),
    legend.text = element_text(size = 8),
    legend.margin = margin(0, 0, 0, 0),
    
    plot.margin = margin(3, 5, 3, 5)
  )

# ---------------------------
# 3. save
# ---------------------------
ggsave(
  file.path(outdir, "13_module_pathway_heatmap_compact.pdf"),
  plot = p_heatmap,
  width = 6.2,
  height = 2.8
)

cat("\n>>> ====== 完成! ======\n")
cat(">>> 输出:\n")
cat("    - 13_module_pathway_heatmap.pdf:     模块-通路热图\n")
cat("    - module_pathway_correlation.csv:    详细结果表\n")

})()

}

writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
message("Original-code workflow
        completed: ", result_dir)
})()

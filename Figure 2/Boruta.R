

# --------------------------- 1. 路径 -----------------------------
input_dir <- "D:/mechanical-related/RA机器学习/RA提取交集基因"
boruta_dir <- "D:/mechanical-related/RA机器学习/Boruta-RA"
output_dir <- file.path(boruta_dir, "corrected_results")

expression_file <- file.path(input_dir, "merged_after_batch_removal.txt")
candidate_file <- file.path(input_dir, "IntersectionGenes.txt")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------- 2. 参数 ------------------------------
expected_gene_number <- 18L
random_seed <- 123L
max_runs <- 100L

# -------------------------- 3. 依赖包 ----------------------------
required_packages <- c("limma", "Boruta", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "缺少以下R包，请先安装：",
    paste(missing_packages, collapse = ", ")
  )
}

# ---------------------- 4. 检查输入文件 --------------------------
if (!file.exists(expression_file)) {
  stop("找不到表达矩阵：", expression_file)
}

if (!file.exists(candidate_file)) {
  stop("找不到候选基因文件：", candidate_file)
}

# -------------------- 5. 读取候选基因列表 ------------------------
# IntersectionGenes.txt当前没有表头，第一行就是CCL2。
# 使用readLines可避免把第一个基因误当作表头。
candidate_lines <- readLines(
  candidate_file,
  warn = FALSE,
  encoding = "UTF-8"
)

# 如文件中出现逗号或制表符，只取每行第一列。
candidate_genes <- vapply(
  strsplit(candidate_lines, "[,\\t]"),
  function(x) x[1],
  character(1)
)

candidate_genes <- sub("^\\ufeff", "", candidate_genes)
candidate_genes <- gsub('^"|"$', "", candidate_genes)
candidate_genes <- trimws(candidate_genes)
candidate_genes <- candidate_genes[
  !is.na(candidate_genes) & nzchar(candidate_genes)
]

# 如果以后文件增加了表头，则自动去掉常见表头名称。
if (
  length(candidate_genes) > 0 &&
  tolower(candidate_genes[1]) %in% c(
    "gene", "genes", "gene_id", "geneid", "symbol", "id"
  )
) {
  candidate_genes <- candidate_genes[-1]
}

duplicated_candidate_genes <- unique(
  candidate_genes[duplicated(candidate_genes)]
)

if (length(duplicated_candidate_genes) > 0) {
  stop(
    "候选基因列表存在重复：",
    paste(duplicated_candidate_genes, collapse = ", ")
  )
}

candidate_genes <- unique(candidate_genes)

cat("候选基因数：", length(candidate_genes), "\n")
cat("候选基因：", paste(candidate_genes, collapse = ", "), "\n")

if (length(candidate_genes) != expected_gene_number) {
  stop(
    "候选基因数应为", expected_gene_number,
    "，实际读取到", length(candidate_genes), "。"
  )
}

if (!"CCL2" %in% candidate_genes) {
  stop("候选基因列表中未检测到CCL2，请检查IntersectionGenes.txt。")
}

# ---------------------- 6. 读取表达矩阵 --------------------------
raw_data <- read.delim(
  expression_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

if (ncol(raw_data) < 2) {
  stop("表达矩阵至少应包含一列基因ID和一列样本表达值。")
}

gene_ids <- trimws(as.character(raw_data[[1]]))
expression_data <- raw_data[, -1, drop = FALSE]

# 明确转换为数值，防止字符型表达值进入分析。
expression_numeric <- lapply(expression_data, function(x) {
  original_missing <- is.na(x) | trimws(as.character(x)) == ""
  converted <- suppressWarnings(as.numeric(as.character(x)))
  newly_missing <- is.na(converted) & !original_missing

  if (any(newly_missing)) {
    stop("表达矩阵中存在无法转换为数值的内容。")
  }

  converted
})

expression_matrix <- as.matrix(
  data.frame(
    expression_numeric,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
)

rownames(expression_matrix) <- gene_ids
colnames(expression_matrix) <- colnames(expression_data)
storage.mode(expression_matrix) <- "double"

if (anyNA(expression_matrix)) {
  stop("表达矩阵中存在缺失值，Boruta运行前需要处理。")
}

# 对重复基因名称取平均值，保留limma::avereps的原有处理逻辑。
expression_averaged <- limma::avereps(
  expression_matrix,
  ID = rownames(expression_matrix)
)

# 保留平均表达量大于0的基因。
expression_filtered <- expression_averaged[
  rowMeans(expression_averaged) > 0,
  ,
  drop = FALSE
]

# ------------------- 7. 严格检查18个基因 -------------------------
missing_genes <- setdiff(candidate_genes, rownames(expression_filtered))

if (length(missing_genes) > 0) {
  stop(
    "以下候选基因未在过滤后的表达矩阵中找到：",
    paste(missing_genes, collapse = ", ")
  )
}

# 按IntersectionGenes.txt中的顺序提取，避免intersect改变解释顺序。
gene_expression <- expression_filtered[
  candidate_genes,
  ,
  drop = FALSE
]

if (nrow(gene_expression) != expected_gene_number) {
  stop(
    "提取后的表达矩阵应包含", expected_gene_number,
    "个基因，实际为", nrow(gene_expression), "。"
  )
}

if (!"CCL2" %in% rownames(gene_expression)) {
  stop("提取后的表达矩阵中缺少CCL2。")
}

# 保存正确的18基因表达矩阵。
gene_expression_output <- data.frame(
  ID = rownames(gene_expression),
  gene_expression,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  gene_expression_output,
  file.path(output_dir, "geneexp_18genes.csv"),
  row.names = FALSE,
  quote = FALSE
)

# ----------------------- 8. 构建分组 -----------------------------
sample_names <- colnames(gene_expression)
group_character <- rep(NA_character_, length(sample_names))

group_character[
  grepl("_Normal$", sample_names, ignore.case = TRUE)
] <- "Normal"

group_character[
  grepl("_RA$", sample_names, ignore.case = TRUE)
] <- "RA"

if (anyNA(group_character)) {
  stop(
    "以下样本名无法根据_Normal或_RA识别分组：",
    paste(sample_names[is.na(group_character)], collapse = ", ")
  )
}

group <- factor(group_character, levels = c("Normal", "RA"))
sample_count <- as.data.frame(table(group), stringsAsFactors = FALSE)
colnames(sample_count) <- c("Group", "Sample_number")

write.csv(
  sample_count,
  file.path(output_dir, "Sample_Counts.csv"),
  row.names = FALSE
)

cat("样本数：\n")
print(sample_count)

# ----------------------- 9. 构建特征矩阵 -------------------------
x_original <- as.data.frame(
  t(gene_expression),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

original_feature_names <- colnames(x_original)
safe_feature_names <- make.names(original_feature_names, unique = TRUE)
feature_name_map <- setNames(original_feature_names, safe_feature_names)
colnames(x_original) <- safe_feature_names

feature_sd <- vapply(x_original, stats::sd, numeric(1), na.rm = TRUE)
zero_variance_features <- names(feature_sd)[
  !is.finite(feature_sd) | feature_sd == 0
]

if (length(zero_variance_features) > 0) {
  stop(
    "以下基因在全部样本中没有变异，无法用于Boruta：",
    paste(feature_name_map[zero_variance_features], collapse = ", ")
  )
}

# 为保持与原分析一致，对各基因进行Z-score标准化。
x_scaled <- as.data.frame(
  scale(x_original),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (anyNA(x_scaled) || any(!is.finite(as.matrix(x_scaled)))) {
  stop("标准化后的特征矩阵包含NA或无穷值。")
}

# ----------------------- 10. 运行Boruta --------------------------
set.seed(random_seed)

boruta_result <- Boruta::Boruta(
  x = x_scaled,
  y = group,
  doTrace = 2,
  maxRuns = max_runs
)

saveRDS(
  boruta_result,
  file.path(output_dir, "Boruta_Result.rds")
)

# ---------------------- 11. 完整状态结果 -------------------------
boruta_stats <- Boruta::attStats(boruta_result)
boruta_stats$Feature_safe_name <- rownames(boruta_stats)
boruta_stats$Gene <- unname(
  feature_name_map[boruta_stats$Feature_safe_name]
)

boruta_stats$decision <- as.character(boruta_stats$decision)
boruta_stats <- boruta_stats[
  order(-boruta_stats$meanImp, -boruta_stats$medianImp),
  ,
  drop = FALSE
]
boruta_stats$Rank_by_meanImp <- seq_len(nrow(boruta_stats))

result_columns <- c(
  "Rank_by_meanImp",
  "Gene",
  "meanImp",
  "medianImp",
  "minImp",
  "maxImp",
  "normHits",
  "decision"
)

boruta_table <- boruta_stats[, result_columns, drop = FALSE]

write.csv(
  boruta_table,
  file.path(output_dir, "Boruta_All_18_Genes_Status_and_Importance.csv"),
  row.names = FALSE
)

status_levels <- c("Confirmed", "Tentative", "Rejected")
status_count <- table(
  factor(boruta_table$decision, levels = status_levels)
)
status_summary <- data.frame(
  Status = status_levels,
  Count = as.integer(status_count),
  stringsAsFactors = FALSE
)

write.csv(
  status_summary,
  file.path(output_dir, "Boruta_Status_Summary.csv"),
  row.names = FALSE
)

# 分别保存三类状态，便于审稿回复核对。
for (status_name in status_levels) {
  status_table <- boruta_table[
    boruta_table$decision == status_name,
    ,
    drop = FALSE
  ]

  write.csv(
    status_table,
    file.path(
      output_dir,
      paste0("Boruta_", status_name, "_Genes.csv")
    ),
    row.names = FALSE
  )
}

# ----------------------- 12. Top 5/Top 10 ------------------------
# Top 5和Top 10只从Confirmed基因中选择。
confirmed_genes <- boruta_table[
  boruta_table$decision == "Confirmed",
  ,
  drop = FALSE
]

if (nrow(confirmed_genes) < 5) {
  stop(
    "Confirmed基因少于5个，不能用Tentative或Rejected基因补足Top 5。"
  )
}

top5_confirmed <- utils::head(confirmed_genes, 5)
top10_confirmed <- utils::head(
  confirmed_genes,
  min(10, nrow(confirmed_genes))
)

write.csv(
  top5_confirmed,
  file.path(output_dir, "Boruta_Top5_Confirmed_Genes.csv"),
  row.names = FALSE
)

write.csv(
  top10_confirmed,
  file.path(output_dir, "Boruta_Top10_Confirmed_Genes.csv"),
  row.names = FALSE
)

# -------------------------- 13. 图形 -----------------------------
pdf(
  file.path(output_dir, "Boruta_Feature_Importance_BasePlot.pdf"),
  width = 12,
  height = 8
)
par(mar = c(9, 4, 4, 2))
plot(
  boruta_result,
  las = 2,
  cex.axis = 0.7,
  xlab = "",
  ylab = "Importance score",
  main = "Boruta feature importance"
)
legend(
  "topleft",
  legend = c(
    "Confirmed",
    "Tentative",
    "Rejected",
    "Shadow features"
  ),
  fill = c("#70AD47", "#FFC000", "#C00000", "#5B9BD5"),
  cex = 0.8,
  bty = "n"
)
dev.off()

plot_data <- boruta_table
plot_data$Gene <- factor(
  plot_data$Gene,
  levels = rev(plot_data$Gene)
)
plot_data$decision <- factor(
  plot_data$decision,
  levels = c("Rejected", "Tentative", "Confirmed")
)

importance_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = Gene, y = meanImp, fill = decision)
) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(
    values = c(
      "Confirmed" = "#70AD47",
      "Tentative" = "#FFC000",
      "Rejected" = "#C00000"
    ),
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "Boruta feature importance for 18 candidate genes",
    subtitle = paste(
      status_summary$Status,
      status_summary$Count,
      collapse = "; "
    ),
    x = NULL,
    y = "Mean importance score (Z-score)",
    fill = "Boruta decision"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "top"
  )

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "Boruta_Feature_Importance_18Genes.pdf"
  ),
  plot = importance_plot,
  width = 8,
  height = 6
)

# -------------------- 14. 参数与软件版本记录 ---------------------
analysis_parameters <- data.frame(
  Parameter = c(
    "Expression_file",
    "Candidate_gene_file",
    "Expected_gene_number",
    "Actual_gene_number",
    "Random_seed",
    "Boruta_maxRuns",
    "Scaling",
    "Ranking_metric",
    "Top5_rule"
  ),
  Value = c(
    expression_file,
    candidate_file,
    expected_gene_number,
    nrow(gene_expression),
    random_seed,
    max_runs,
    "Z-score by gene across samples",
    "meanImp descending among Confirmed genes",
    "Five highest-ranked Confirmed genes"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  analysis_parameters,
  file.path(output_dir, "Boruta_Analysis_Parameters.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "R_SessionInfo.txt")
)

summary_lines <- c(
  "RA Boruta corrected analysis summary",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Candidate genes: ", length(candidate_genes)),
  paste0("Extracted genes: ", nrow(gene_expression)),
  paste0("CCL2 included: ", "CCL2" %in% rownames(gene_expression)),
  paste0("Normal samples: ", sum(group == "Normal")),
  paste0("RA samples: ", sum(group == "RA")),
  paste0(
    "Confirmed/Tentative/Rejected: ",
    paste(status_summary$Count, collapse = "/")
  ),
  paste0("Top 5 confirmed genes: ", paste(top5_confirmed$Gene, collapse = ", ")),
  "Top 5 is a fixed presentation/prioritization rule, not a statistical significance threshold."
)

writeLines(
  summary_lines,
  con = file.path(output_dir, "Boruta_Analysis_Summary.txt")
)

# ----------------------- 15. 控制台汇总 --------------------------
cat("\nBoruta分析完成。\n")
cat("结果目录：", output_dir, "\n")
cat("\n状态统计：\n")
print(status_summary)
cat("\nTop 5 Confirmed基因：\n")
print(top5_confirmed)


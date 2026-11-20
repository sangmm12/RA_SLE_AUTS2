# ==== ROC 分析：联合模型与单基因 ====
# 基于原始 pipeline 精简；保留 rms::lrm 模型和原 ROC 绘图样式。
# 输入：表达矩阵为基因行、样本列，首列是基因名；基因列表无表头。
# 样本名须以 _0 / _1 结尾，0 = 对照，1 = 疾病；请核实实际编码。
# 联合模型在当前数据拟合并评价，输出为训练集（表观）ROC。
# 单基因方向处理沿用原代码，AUC 不代表独立验证结果。

# ============================ 加载必要包 =============================
# 首次使用若未安装：install.packages(c("rms", "pROC"))
library(rms)
library(pROC)

# Step 1: 输入与输出目录
working_directory <- "D:\\mechanical-related\\6RA\\2\\RA"
if (!dir.exists(working_directory)) {
  stop(paste("指定的工作目录不存在：", working_directory))
}
input_expr <- file.path(working_directory, "Sample Type Matrix.csv")
gene_file <- file.path(working_directory, "IntersectionGenes.txt")
# 沿用原目录及 ROC 文件名，重复运行会覆盖同名输出。
result_dir <- file.path(working_directory, "Analysis_Result")
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

# 用函数隔离工作目录及 datadist 设置，结束或报错时恢复原设置。
run_roc_analysis <- function() {
  old_wd <- getwd()
  old_options <- options("datadist")
  on.exit({
    options(old_options)
    setwd(old_wd)
  }, add = TRUE)
  setwd(result_dir)
  cat("==== ROC 分析启动 ====\n")

  # Step 2: 读取表达矩阵和基因列表
  if (!file.exists(input_expr)) stop(sprintf("未找到文件: %s", input_expr))
  if (!file.exists(gene_file)) stop(sprintf("未找到文件: %s", gene_file))
  dat_expr <- read.table(
    input_expr, sep = ",", header = TRUE, row.names = 1,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  gene_tbl <- read.table(
    gene_file, header = FALSE, check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(gene_tbl) < 1) stop("基因文件无有效列")
  feature_genes <- trimws(as.character(gene_tbl[, 1]))
  feature_genes <- unique(feature_genes[!is.na(feature_genes) & feature_genes != ""])
  if (length(feature_genes) < 2) stop("特征基因数过少，联合分析至少需要两个基因")

  # Step 3: 统一基因名并筛选
  # 修正原脚本：修改行名，而不是按替换后的名称重新索引表达矩阵。
  normalized_genes <- gsub("-", "_", rownames(dat_expr))
  if (anyDuplicated(normalized_genes)) {
    stop("统一基因名后存在重复行名，请先确认并处理重复基因")
  }
  rownames(dat_expr) <- normalized_genes
  feature_genes <- unique(gsub("-", "_", feature_genes))
  found_genes <- feature_genes %in% rownames(dat_expr)
  if (!any(found_genes)) stop("IntersectionGenes.txt 中的基因在表达矩阵中均找不到")
  if (any(!found_genes)) {
    warning(paste0("以下基因不在表达矩阵中：",
                   paste(feature_genes[!found_genes], collapse = ",")))
  }
  dat_expr_filt <- dat_expr[feature_genes[found_genes], , drop = FALSE]
  if (nrow(dat_expr_filt) < 2) stop("实际匹配到的基因不足两个")
  if (!all(vapply(dat_expr_filt, is.numeric, logical(1)))) {
    stop("表达矩阵包含非数值列，请检查输入文件")
  }
  if (any(!is.finite(as.matrix(dat_expr_filt)))) {
    stop("所选基因的表达数据含 NA/NaN/Inf，请先处理缺失或非有限数值")
  }
  if (any(vapply(seq_len(nrow(dat_expr_filt)), function(i) {
    length(unique(unlist(dat_expr_filt[i, ], use.names = FALSE))) < 2
  }, logical(1)))) {
    stop("所选基因中有表达值恒定的基因，请先检查")
  }
  cat(sprintf("实际纳入 %d 个基因\n", nrow(dat_expr_filt)))

  # Step 4: 转置并提取分组
  df_expr <- as.data.frame(t(dat_expr_filt))
  sample_names <- rownames(df_expr)
  if (!all(grepl("_[01]$", sample_names))) {
    stop("样本名必须以 _0 或 _1 结尾；请核实 0=对照、1=疾病")
  }
  groups <- sub("^.*_", "", sample_names)
  if (!setequal(unique(groups), c("0", "1"))) {
    stop("必须同时包含 0 组和 1 组样本")
  }
  df_expr$GroupType <- as.integer(groups)
  model_vars <- setdiff(colnames(df_expr), "GroupType")
  if ("GroupType" %in% rownames(dat_expr_filt)) stop("基因名不能与分组列 GroupType 重名")
  if (any(grepl("`", model_vars, fixed = TRUE))) stop("基因名中含反引号，请先检查")
  cat(sprintf("对照组：%d；疾病组：%d\n",
              sum(df_expr$GroupType == 0), sum(df_expr$GroupType == 1)))

  # Step 5: 保留原 lrm 联合模型，生成 ROC 所需预测概率
  ddinfo <- datadist(df_expr)
  # 给 rms 可访问的数据分布对象，避免依赖全局工作区中的旧对象。
  options(datadist = ddinfo)
  reg_formula <- as.formula(paste(
    "GroupType ~", paste(paste0("`", model_vars, "`"), collapse = " + ")
  ))
  lrm_fit <- lrm(reg_formula, data = df_expr, x = TRUE, y = TRUE)
  if (isTRUE(lrm_fit$fail)) stop("lrm 模型拟合失败，请检查共线性、样本量或完全分离")
  pred_probs <- as.numeric(predict(lrm_fit, newdata = df_expr, type = "fitted"))
  if (length(pred_probs) != nrow(df_expr) || any(!is.finite(pred_probs))) {
    stop("预测概率无效或与样本数不一致")
  }

  # 保留用于追溯 ROC 的样本预测表。
  df_out <- data.frame(
    Sample = rownames(df_expr), GroupType = df_expr$GroupType,
    Predicted_Prob = pred_probs
  )
  write.csv(df_out, "Sample_Predicted_Probabilities.csv", row.names = FALSE)

  # Step 6: 绘制 ROC 并输出 AUC（绘图参数沿用原代码）
  # ==== 联合模型 ROC（主模型红色，风格与单基因统一）====


  roc_y <- as.numeric(df_expr$GroupType)
  roc_pred <- pred_probs
  roc_obj <- roc(roc_y, roc_pred, levels = c(0,1), direction = "<")

  # 主模型/单基因统一色板（前几个色都与后面单基因一致）
  my_cols <- c("red", "deepskyblue", "forestgreen", "orange", 
               "purple", "gray40", "black", "magenta", "gold", "brown")
  while(length(my_cols) < length(model_vars) + 1) { 
    my_cols <- c(my_cols, rainbow(length(model_vars) + 1 - length(my_cols)))
  }

  ## ====== 联合主模型ROC ======
  pdf("CombinedGenes_ROC_SCI.pdf", width = 6, height = 6)
  par(mar = c(5,6,4,2)+0.1, cex = 1.3)
  plot(1-roc_obj$specificities, roc_obj$sensitivities, type = "l",
       col = my_cols[1], lwd = 4, lty = 1,
       xlab = expression("1 - Specificity"), ylab = "Sensitivity",
       main = " Model ROC Curve", 
       cex.lab = 1.4, cex.axis = 1.15, cex.main = 1.45,
       xlim = c(0,1), ylim = c(0,1))
  abline(0, 1, lty = 2, col = "gray70", lwd = 2)
  auc_val <- as.numeric(auc(roc_obj))
  legend("bottomright", legend = sprintf("AUC = %.3f", auc_val), col = my_cols[1], 
         lwd = 4, lty=1, bty = "n", cex = 1.2)
  dev.off()
  cat(">> 合并主模型 ROC 曲线 (红色) 已输出至 CombinedGenes_ROC_SCI.pdf\n")
  # ==== 单基因多色ROC优化版 ====
  roc_list <- list(); auc_list <- c(); leglab <- c(); auc_dir <- c()


  pdf("IndividualGenes_ROC_SCI.pdf", width = 8, height = 8)
  par(mar = c(5, 6, 4, 2)+0.1, cex = 1.3)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
       xlab = expression("1 - Specificity"), ylab = "Sensitivity",
       main = "Gene ROC Curves", 
       cex.lab = 1.4, cex.axis = 1.15, cex.main = 1.45
  )
  abline(0, 1, lty = 2, col = "gray70", lwd = 2)

  for(i in seq_along(model_vars)) {
    g <- model_vars[i]
    # 先用direction="auto"自动判断
    cur_roc <- roc(roc_y, df_expr[[g]], levels = c(0,1), direction = "auto")
    cur_auc <- as.numeric(auc(cur_roc))
  
    # 沿用原脚本的方向反转规则；外部验证时应固定训练集确定的方向。
    flip <- FALSE
    if(cur_auc < 0.5) {
      flip <- TRUE
      cur_auc <- 1 - cur_auc
      # 必须重新画ROC曲线，用反向direction
      cur_roc <- roc(roc_y, df_expr[[g]], levels = c(0,1), direction = ifelse(cur_roc$direction==">", "<", ">"))
    }
  
    lines(1-cur_roc$specificities, cur_roc$sensitivities, col = my_cols[i+1], lwd = 3)
    roc_list[[g]] <- cur_roc
    auc_list <- c(auc_list, cur_auc)
  
    # 写入方向注释
    mean0 <- mean(df_expr[[g]][roc_y==0], na.rm=TRUE)
    mean1 <- mean(df_expr[[g]][roc_y==1], na.rm=TRUE)
    direction_note <- ifelse(mean1 > mean0, paste0("(组1高)"), paste0("(组0高)"))
    auc_dir <- c(auc_dir, direction_note)
  
    # 图例写明方向和AUC
    leglab <- c(leglab, sprintf("%s  AUC=%.3f", g, cur_auc))
  
  }
  legend("bottomright", legend = leglab, col = my_cols[2:(length(model_vars)+1)], 
         lwd = 3, lty = 1, bty = "n", cex = 1.15, y.intersp = 1.18)
  dev.off()

  # AUC同步输出加强版
  gene_auc_df <- data.frame(
    Gene = model_vars,
    AUC = auc_list,
    Direction = auc_dir, ROC_Direction = vapply(roc_list, function(x) x$direction, character(1))
  )
  write.csv(gene_auc_df, file = "IndividualGenes_AUC.csv", row.names = FALSE)
  cat(">> 单基因 ROC 曲线和 AUC 已输出至 IndividualGenes_ROC_SCI.pdf / IndividualGenes_AUC.csv\n")


  # 联合模型 AUC 汇总
  write.csv(
    data.frame(Model = "CombinedGenes", AUC = auc_val),
    "CombinedGenes_AUC.csv", row.names = FALSE
  )
  cat("==== ROC 分析完成 ====\n")
  cat("输出目录：", getwd(), "\n")
  invisible(list(model = lrm_fit, combined_roc = roc_obj, gene_rocs = roc_list))
}

roc_results <- run_roc_analysis()

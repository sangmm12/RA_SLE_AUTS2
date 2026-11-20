# --- 加载必要包 ---
library(ggvenn)
library(grid)
library(gridExtra)
library(ggplot2)
library(eulerr)

# --- 设置工作目录 ---
workDir <- "D:\\mechanical-related\\KOOEVEEN\\RA-OE ∩ SLE-OE"
setwd(workDir)

# --- 创建输出文件夹 ---
outputDir <- "output_folder"
if (!dir.exists(outputDir)) dir.create(outputDir)

# --- 自动识别目录下的CSV和TXT文件 ---
files_txt <- list.files(pattern = "\\.txt$", ignore.case = TRUE)
files_csv <- list.files(pattern = "\\.csv$", ignore.case = TRUE)

# 排除输出文件
exclude_files <- c("compound.txt", "IntersectionGenes.txt")
files_txt <- files_txt[!files_txt %in% exclude_files]

files <- c(files_txt, files_csv)
cat(sprintf("发现 %d 个TXT文件，%d 个CSV文件\n", length(files_txt), length(files_csv)))
cat("文件列表：\n")
print(files)

# --- 读取所有文件中的基因信息，整理到geneList ---
geneList <- list()
for (inputFile in files) {
  cat(sprintf("\n处理文件: %s\n", inputFile))

  # 判断文件类型
  file_ext <- tolower(tools::file_ext(inputFile))

  if (file_ext == "csv") {
    # CSV文件：有表头，使用逗号分隔
    rt <- read.csv(inputFile, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    geneNames <- as.character(rt[, 1])
    cat(sprintf("  CSV文件，有表头，读取第一列\n"))
  } else {
    # TXT文件：无表头，使用制表符或空格分隔
    rt <- read.table(inputFile, header = FALSE, stringsAsFactors = FALSE,
                     check.names = FALSE, fill = TRUE)
    geneNames <- as.character(rt[, 1])
    cat(sprintf("  TXT文件，无表头，读取第一列\n"))
  }

  # 清理基因名
  geneNames <- trimws(geneNames)           # 去掉前后空格
  geneNames <- geneNames[geneNames != ""]  # 去掉空字符串
  geneNames <- geneNames[!is.na(geneNames)] # 去掉NA
  uniqGene <- unique(geneNames)

  # 使用文件名（去扩展名）作为集合名
  setName <- tools::file_path_sans_ext(basename(inputFile))
  geneList[[setName]] <- uniqGene
  cat(sprintf("  共提取 %d 个不重复基因\n", length(uniqGene)))
}

# --- 显示各集合的基因数量 ---
cat("\n各集合基因数量：\n")
for (name in names(geneList)) {
  cat(sprintf("  %s: %d 个基因\n", name, length(geneList[[name]])))
}

# --- 计算所有基因名的交集 ---
intersectionGenes <- Reduce(intersect, geneList)
intersectionCount <- length(intersectionGenes)
cat(sprintf("\n交集基因数量: %d\n", intersectionCount))
cat("交集基因列表:\n")
print(intersectionGenes)

# --- 绘制Venn图并保存 ---
# 定义颜色（支持最多5个集合）
colors <- c("#FFA700", "#1E90FF", "#4DAF4A", "#984EA3", "#FF7F00")
n_sets <- length(geneList)

pdf(file = file.path(outputDir, "venn.pdf"), width = 11, height = 8)
venn_plot <- ggvenn(
  geneList,
  show_percentage = TRUE,
  stroke_color = "white",
  stroke_size = 0.5,
  fill_color = colors[1:n_sets],
  set_name_color = colors[1:n_sets],
  set_name_size = 7,
  text_size = 5
)
# 添加交集个数
grid.arrange(venn_plot,
             textGrob(paste("Intersection Gene Count:", intersectionCount),
                      gp = gpar(fontsize = 14, fontface = "bold")),
             ncol = 1, heights = c(5, 1))
dev.off()
cat("Venn图已保存到 output_folder/venn.pdf\n")

# --- 输出交集基因到output_folder ---
write.table(
  intersectionGenes,
  file = file.path(outputDir, "IntersectionGenes.txt"),
  sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE
)
cat("交集基因已保存到 output_folder/IntersectionGenes.txt\n")

# --- 绘制各集合基因数量柱状图 ---
setSizes <- sapply(geneList, length)
pdf(file = file.path(outputDir, "set_barplot.pdf"), width = 8, height = 5)
ggplot(data.frame(Set = names(setSizes), Size = setSizes), aes(x = Set, y = Size, fill = Set)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Size), vjust = -0.5, size = 4) +
  scale_fill_manual(values = colors[1:n_sets]) +
  theme_minimal(base_size = 14) +
  ylab("Gene Count") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  ggtitle("Gene Count per Set")
dev.off()
cat("柱状图已保存到 output_folder/set_barplot.pdf\n")

# --- 绘制Euler图 ---
fit <- euler(geneList)
pdf(file = file.path(outputDir, "euler.pdf"), width = 8, height = 8)
plot(fit,
     fills = list(fill = colors[1:n_sets], alpha = 0.7),
     labels = list(font = 2, cex = 1.2),
     edges = list(lwd = 2),
     quantities = list(cex = 1.1))
dev.off()
cat("Euler图已保存到 output_folder/euler.pdf\n")

# --- 完成提示 ---
cat("\n========================================\n")
cat("处理完成！\n")
cat("========================================\n")
cat(sprintf("输入文件数: %d\n", length(files)))
cat(sprintf("交集基因数: %d\n", intersectionCount))
cat("输出文件:\n")
cat("  - output_folder/venn.pdf: Venn图\n")
cat("  - output_folder/euler.pdf: Euler图\n")
cat("  - output_folder/set_barplot.pdf: 各集合基因数量柱状图\n")
cat("  - output_folder/IntersectionGenes.txt: 交集基因列表\n")
cat("========================================\n")


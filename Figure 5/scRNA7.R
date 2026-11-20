
#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(ROGUE)
library(clustree)
library(harmony)
library(SingleR)
#SingleR version 2.4.1

#目录
setwd("D:/mechanical-related/RA-scRNA/9")

#读取
pbmc = readRDS("D:/mechanical-related/RA-scRNA/9/5pbmc_UMPA.TSNE.rds")

#加载参考数据库
load("ref_Human_all.RData")

#获取表达矩阵
testdata = GetAssayData(object = pbmc@assays$RNA, layer = "counts")
#testdata = GetAssayData(object = pbmc@assays$SCT, layer = "data")

#获取clusters
clusters <- pbmc@meta.data$seurat_clusters

#参考数据库的大类和小类细胞群
table(ref_Human_all@colData@listData[["label.main"]])
table(ref_Human_all@colData@listData[["label.fine"]])

#运行singleR
cellpred <- SingleR(test = testdata, ref = ref_Human_all, clusters = clusters, assay.type.test = "logcounts", 
                    labels = ref_Human_all@colData@listData[["label.main"]], assay.type.ref = "logcounts")

#获取每个cluster的细胞类型
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

#添加到seurat.metadata对象中
pbmc@meta.data$SingleR = "NA"
for(i in 1:nrow(celltype)){
  pbmc@meta.data[which(pbmc$seurat_clusters == celltype$ClusterID[i]),'SingleR'] <- celltype$celltype[i]
}

#画图
p = plotScoreHeatmap(cellpred)
ggsave("1.pdf", p, width = 12, height = 5)

p1 <- DimPlot(pbmc, group.by = "SingleR", label = T,reduction = "tsne")
p2 <- DimPlot(pbmc, group.by = "SingleR", label = T,reduction = "umap")
p <- p1 | p2
ggsave("2.pdf", p, width = 12, height = 5)

#保存
saveRDS(pbmc,"6pbmc_SingleR.rds")

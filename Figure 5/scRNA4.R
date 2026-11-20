
#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)

#目录
setwd("D:/mechanical-related/RA-scRNA/6")

#读取
pbmc = readRDS("D:/mechanical-related/RA-scRNA/6/2pbmc_double.rds")

#细胞周期评分
pbmc <- NormalizeData(pbmc)

#获取G2M期相关基因
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(pbmc))

#获取S期相关基因
s_genes <- cc.genes$s.genes    
s_genes <- CaseMatch(search=s_genes, match=rownames(pbmc))

#细胞周期阶段评分
pbmc <- CellCycleScoring(pbmc, g2m.features=g2m_genes, s.features=s_genes)

colnames(pbmc@meta.data)
table(pbmc$Phase)

#画图
DimPlot(pbmc,group.by = "Phase",reduction = "tsne")

#保存
saveRDS(pbmc,"3pbmc_CellCycle.rds")



#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)

#目录
setwd("D:/mechanical-related/RA-scRNA/5")

#读取
pbmc = readRDS("D:/mechanical-related/RA-scRNA/5/1pbmc_qc.rds")

#%>%（向右操作符，forward-pipe operator）是最常用的一种操作符，
#就是把左侧准备的数据或表达式，传递给右侧的函数调用或表达式进行运行，
#可以连续操作就像一个链条一样。

#数据标准化、寻找可变特征、数据缩放
pbmc = pbmc %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

#PCA、UMAP、寻找邻居、聚类
pbmc = pbmc %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>%  
  RunTSNE(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.1)

#首先获得最佳的pK值
#pK表示领域大小
sweep.res.list <- paramSweep(pbmc, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_best = bcmvn %>% 
  dplyr::arrange(desc(BCmetric)) %>% 
  dplyr::pull(pK) %>% 
  .[1] %>% as.character() %>% as.numeric()

#然后估算出双细胞群中，homotypic doublets的比例
annotations <- pbmc$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
print(homotypic.prop)

#双细胞占比为7%左右
nExp_poi <- round(0.07*nrow(pbmc@meta.data))        
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 

# 模拟出的artificial doublet数量。不同取值对识别结果影响不大，默认为0.25
pbmc <- doubletFinder(pbmc, PCs = 1:30, 
                       pN = 0.25, pK = pk_best, nExp = nExp_poi.adj, 
                       reuse.pANN = NULL, sct = FALSE)

#若出现：错误与xtfrm.data.frame(x)：无法xtfrm数据帧
#将上述的reuse.pANN＝FALSE改成NULL

#将列名改为"Double_score"和"Is_Double"
colnames(pbmc@meta.data)

colnames(pbmc@meta.data)[length(colnames(pbmc@meta.data))-1] <- "Double_score"
colnames(pbmc@meta.data)[length(colnames(pbmc@meta.data))] <- "Is_Double"

#查看DoubletFinder分析结果
head(pbmc@meta.data[, c("Double_score", "Is_Double")])

#绘制DoubletFinder分类的tsne图
DimPlot(pbmc, reduction = "tsne", group.by = "Is_Double")

#绘制双细胞分类的小提琴图
VlnPlot(pbmc, group.by = "Is_Double", 
        features = c("nCount_RNA", "nFeature_RNA"), 
        pt.size = 0, ncol = 2)

#过滤非单细胞数据
#pbmc <- subset(pbmc, Is_Double == "Singlet")

#保存
saveRDS(pbmc,"2pbmc_double.rds")




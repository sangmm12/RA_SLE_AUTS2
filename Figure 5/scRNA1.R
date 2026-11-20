
#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)

#工作目录
setwd("D:/mechanical-related/RA-scRNA/3")


#####1.matrix.mtx、genes.tsv和barcodes.tsv####

#目录
setwd("D:/mechanical-related/RA-scRNA/3/GSE200815")

# 获取数据文件夹下的所有样本文件列表
samples <- list.files("seurat/")
samples

# 创建一个空的列表
seurat_list <- list()

#读取数据并创建Seurat对象
#删除，小于200个基因表达的细胞，小于3个细胞表达的基因
for (sample in samples) {
  #文件路径
  data.path <- paste0("seurat/", sample)
  
  #读取10x数据
  seurat_data <- Read10X(data.dir = data.path)

  #创建Seurat对象
  seurat_obj <- CreateSeuratObject(counts = seurat_data,project = sample,min.features = 200,min.cells = 3)
  
  #添加到列表中
  seurat_list <- append(seurat_list, seurat_obj)
}

#合并
seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)

#将Layers融合
pbmc = JoinLayers(seurat_combined)

saveRDS(pbmc, file = "GSE200815_merged_seurat.rds")


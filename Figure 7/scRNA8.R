
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
library(dplyr)

#目录
setwd("D:/mechanical-related/7-SLE-scRNA/10")

#读取
pbmc = readRDS("D:/mechanical-related/7-SLE-scRNA/10/6pbmc_SingleR.rds")

DefaultAssay(pbmc)
#DefaultAssay(pbmc) = "RNA"
table(Idents(pbmc))
#Idents(pbmc) = "nCount_RNA"

#寻找marker
pbmc.markers1 <- FindAllMarkers(pbmc, only.pos = TRUE,logfc.threshold = 1,)
write.csv(pbmc.markers1,file="markers.1.SCT.csv")

DefaultAssay(pbmc) = "RNA"
pbmc.markers2 <- FindAllMarkers(pbmc, only.pos = TRUE,logfc.threshold = 1)
write.csv(pbmc.markers2,file="markers.2.RNA.csv")
DefaultAssay(pbmc)
DefaultAssay(pbmc) = "RNA"

#创建marker集合
markers <- c("PTPRC", #immune
             "KRT19", #epithelial
             "COL1A1","MME","DCN","PECAM1") #stromal

#画图
p <- FeaturePlot(pbmc, features = markers, ncol = 3)
ggsave("1.pdf", p, width = 15, height = 10)

p <- DotPlot(pbmc, features = markers) + RotatedAxis()
ggsave("2.pdf", p, width = 14, height = 6)

p <- VlnPlot(pbmc, features = markers, stack = T, flip = T) + NoLegend()
ggsave("3.pdf", p, width = 14, height = 6)

pdf(file="4.pdf",width=7,height=6)
DimPlot(pbmc, reduction = "umap", label = T)
dev.off()

#第一次
#immune epithelial stromal
pbmc$celltype.1 <- recode(pbmc@meta.data$seurat_clusters,
                             "0" = "immune",
                             "1" = "immune",
                             "2" = "immune",
                             "3" = "immune",
                             "4" = "immune",
                             "5" = "immune",
                             "6" = "immune",
                             "7" = "immune",
                             "8" = "immune",
                             "9" = "immune",
                             "10" = "immune",
                             "11" = "immune",
                             "12" = "immune",
                             "13" = "immune",
                             "14" = "immune",
                             "15" = "immune",
                             "16" = "immune",
                             "17" = "immune",
                             "18" = "immune",
                             "19" = "immune",
                             "20" = "immune",
                             "21" = "immune",
                             "22" = "immune",
                             "23" = "immune")

#查看各个细胞群的亚型
table(pbmc@meta.data$celltype.1)

#换个配色
Biocols = c('#AB3282', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
            '#E95C59', '#E59CC4','#E5D2DD' , '#23452F', '#BD956A', '#8C549C', '#585658',
            '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
            '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
            '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
            '#968175')

#画图
p = DimPlot(pbmc, reduction = "umap", label = T,group.by = "celltype.1",cols = Biocols)
ggsave("5.pdf", p, width = 7, height = 6)







#第二次
#Marker基因标记
markers <- c("CD3E","IL7R","CCR7", #T_cell
             "LYZ","CD14","CD68","S100A9","FCGR3A", #monocyte
             "KLRB1","NKG7","GNLY", #NK_cell
             "CD19","CD79A","MS4A1", #B_cell
             "CD1C","FCER1A","CST3","LILRA4","CLEC4C") #DC_cell

#画图
# p <- FeaturePlot(pbmc, features = markers, ncol = 3)
# ggsave("6.pdf", p, width = 15, height = 32)
p <- DotPlot(pbmc, features = markers) + RotatedAxis()
ggsave("7.pdf", p, width = 14, height = 6)
p <- VlnPlot(pbmc, features = markers, stack = T, flip = T) + NoLegend()
ggsave("8.pdf", p, width = 14, height = 6)

#Fibroblast  Endothelia  Mast  Luminal Basal/intermediate  Monolytic  T
pbmc$celltype.main <- recode(pbmc@meta.data$seurat_clusters,
                             "0" = "T_cell",
                             "1" = "T_cell",
                             "2" = "monocyte",
                             "3" = "T_cell",
                             "4" = "T_cell",
                             "5" = "NK_cell",
                             "6" = "T_cell",
                             "7" = "NK_cell",
                             "8" = "NK_cell",
                             "9" = "monocyte",
                             "10" = "B_cell",
                             "11" = "monocyte",
                             "12" = "NK_cell",
                             "13" = "T_cell",
                             "14" = "T_cell",
                             "15" = "T_cell",
                             "16" = "B_cell",
                             "17" = "NK_cell",
                             "18" = "monocyte",
                             "19" = "monocyte",
                             "20" = "B_cell",
                             "21" = "DC_cell",
                             "22" = "DC_cell",
                             "23" = "monocyte")
#查看各个细胞群的亚型
table(pbmc@meta.data$celltype.main)
#画图
p = DimPlot(pbmc, reduction = "umap", label = T,group.by = "celltype.main")
ggsave("9.pdf", p, width = 7, height = 6)

p1 = DimPlot(pbmc, reduction = "umap", label = T,group.by = "celltype.main")
p2 = DimPlot(pbmc, reduction = "umap", label = T,group.by = "celltype.1")
ggsave("10.pdf", p1|p2, width = 15, height = 6)









#第四次
#Marker基因标记
markers <- c("CD3E","IL7R","CCR7", #T_cell
             "LYZ","CD14","CD68","S100A9", #monocyte
             "KLRB1","NKG7","GNLY", #NK_cell
             "CD19","CD79A","MS4A1", #B_cell
             "CD1C","FCER1A","CST3","LILRA4") #DC_cell
             

#画图
p <- DotPlot(pbmc, features = markers,group.by = "celltype.main") + RotatedAxis()
ggsave("14.pdf", p, width = 14, height = 6)

#查看各个细胞群的亚型
table(pbmc@meta.data$celltype.main)









#分组
table(Idents(pbmc))
Idents(pbmc) = "celltype.main"
table(Idents(pbmc))

#矩阵
DefaultAssay(pbmc)

#寻找每个细胞类型的markers
pbmc.markers3 <- FindAllMarkers(pbmc, only.pos = TRUE,logfc.threshold = 1,min.pct = 0.3)
write.csv(pbmc.markers3,file="markers.celltype.RNA.csv")

#提取前10的marker
top10 <- pbmc.markers3 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#对这些基因进行scaledata
markers = as.data.frame(top10[,"gene"])
pbmc <- ScaleData(pbmc, features = as.character(unique(markers$gene)))

#画图
p = DoHeatmap(pbmc,
          features = as.character(unique(markers$gene)),
          group.by = "celltype.main")
ggsave("15.pdf", p, width = 10, height = 9)

#每个细胞亚群抽最少的细胞类型
allCells = names(Idents(pbmc))
allType = levels(Idents(pbmc))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(pbmc)== x ]
  cg=sample(cgCells,min(table(pbmc@meta.data$celltype.main)))
  cg
}))

#提取
cg_sce = pbmc[, allCells %in% choose_Cells]
table(Idents(cg_sce))

#画图
p = DoHeatmap(cg_sce,
              features = as.character(unique(markers$gene)),
              group.by = "celltype.main")
ggsave("16.pdf", p, width = 10, height = 9)




#堆叠柱状图
cell.prop<-as.data.frame(prop.table(table(pbmc@meta.data$celltype.main, pbmc@meta.data$orig.ident)))
colnames(cell.prop)<-c("cluster","group","proportion")

p = ggplot(cell.prop,aes(group,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))
ggsave("17.pdf", p, width = 10, height = 9)

#保存
saveRDS(pbmc,"7pbmc_celltype.rds")

#提取细胞类型进行细分亚群定义
pbmc.T = pbmc[, Idents(pbmc) %in% c("T")]


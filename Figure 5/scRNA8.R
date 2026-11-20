
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
setwd("D:/mechanical-related/RA-scRNA/10")

#读取
pbmc = readRDS("D:/mechanical-related/RA-scRNA/10/6pbmc_SingleR.rds")

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
                             "0" = "stromal",
                             "1" = "immune",
                             "2" = "stromal",
                             "3" = "stromal",
                             "4" = "stromal",
                             "5" = "stromal",
                             "6" = "stromal",
                             "7" = "stromal",
                             "8" = "immune",
                             "9" = "stromal",
                             "10" = "stromal",
                             "11" = "immune",
                             "12" = "stromal",
                             "13" = "stromal",
                             "14" = "stromal",
                             "15" = "stromal",
                             "16" = "stromal",
                             "17" = "immune",
                             "18" = "immune",
                             "19" = "immune",
                             "20" = "stromal",
                             "21" = "stromal",
                             "22" = "immune",
                             "23" = "stromal",
                             "24" = "immune",
                             "25" = "immune",
                             "26" = "stromal",
                             "27" = "stromal")

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
markers <- c("MMP3","PLA2G2A","DCN","COL3A1","PRG4", #Fibroblast
             "SELE","TM4SF1","VWF","ACKR1","PECAM1", #Endothelial
             "IGHG1","MZB1","SDC1", #Plasma
             "C1QA","IL1B","C1QB","RNASE1","LYZ", #macrophage
             "CD19","CD79A","MS4A1", #B_cell
             "IL7R","CD3D","CD3E","CD8A", #T_cell
             "FGFBP2","FCGR3A","CX3CR1") #NK_cell

#画图
# p <- FeaturePlot(pbmc, features = markers, ncol = 3)
# ggsave("6.pdf", p, width = 15, height = 32)
p <- DotPlot(pbmc, features = markers) + RotatedAxis()
ggsave("7.pdf", p, width = 14, height = 6)
p <- VlnPlot(pbmc, features = markers, stack = T, flip = T) + NoLegend()
ggsave("8.pdf", p, width = 14, height = 6)

#Fibroblast  Endothelia  Mast  Luminal Basal/intermediate  Monolytic  T
pbmc$celltype.main <- recode(pbmc@meta.data$seurat_clusters,
                             "0" = "Fibroblast",
                             "1" = "T_cell",
                             "2" = "Fibroblast",
                             "3" = "Endothelial",
                             "4" = "Fibroblast",
                             "5" = "Endothelial",
                             "6" = "Fibroblast",
                             "7" = "Fibroblast",
                             "8" = "T_cell",
                             "9" = "Endothelial",
                             "10" = "Fibroblast",
                             "11" = "macrophage",
                             "12" = "Fibroblast",
                             "13" = "Endothelial",
                             "14" = "Endothelial",
                             "15" = "Fibroblast",
                             "16" = "Endothelial",
                             "17" = "macrophage",
                             "18" = "macrophage",
                             "19" = "T_cell",
                             "20" = "Endothelial",
                             "21" = "Fibroblast",
                             "22" = "Plasma",
                             "23" = "Fibroblast",
                             "24" = "Plasma",
                             "25" = "B_cell",
                             "26" = "Endothelial",
                             "27" = "Endothelial")
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
markers <- c("MMP3","PLA2G2A","DCN","COL3A1","PRG4", #Fibroblast
             "SELE","TM4SF1","VWF","ACKR1","PECAM1", #Endothelial
             "IGHG1","MZB1","SDC1", #Plasma
             "C1QA","IL1B","C1QB","RNASE1","LYZ", #macrophage
             "CD79A","MS4A1", #B_cell
             "IL7R","CD3D","CD3E","CD8A") #T_cell
             

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




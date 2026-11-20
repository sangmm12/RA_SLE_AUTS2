library(CellChat)

rm(list=ls())
library(tidyverse)
library(Seurat)
setwd()
load() #此处加载完成细胞注释的数据

Idents(scRNA) #查看默认分组
Idents(scRNA) <- scRNA@meta.data$celltype.main #默认分组改为细胞类型
Idents(scRNA) #查看默认分组
table(scRNA$group) #查看group分组情况

#

#1.创建cellchat对象####

#1.方法1:直接创建cellchat
cellchat <- createCellChat(object   = scRNA,     # seurat格式文件名
                           group.by = "ident",   # 默认Idents(scRNA)的当前分组
                           assay    = "RNA"      # 使用RNA assay数据
)



#选择CellChatDB:CellChatDB.human人 CellChatDB.mouse鼠
CellChatDB <- CellChatDB.human #人
#显示数据库的结构
showDatabaseCategory(CellChatDB)

#2.数据库的选择####
#1.只使用分泌型信号
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")


#最终选定数据库输入cellchat
cellchat@DB <- CellChatDB.use

#3.数据预处理及多线程设置####
#数据过滤
cellchat <- subsetData(cellchat) 
#多线程设置
future::plan("multisession", workers = 6) #根据自己的CPU线程数选择合适线程;multisession/multicore
#鉴定过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
#鉴定过表互作
cellchat <- identifyOverExpressedInteractions(cellchat)

#如果互作数量少(如number<150),可进一步将基因表达数据投影到PPI(蛋白互作网络)
#直接分析：只看配体IL6+受体IL6R的表达,如果其中一个没表达,就判定这条通路无互作；
#投影到PPI后:还会看IL6的一些互作蛋白的表达,只要这些互作蛋白有表达就认为IL6通路有互作。
#注意：运行PPI时,应在函数computeCommunoProb中设置raw.use=FALSE,以便使用投影数据
# cellchat <- projectData(cellchat, PPI.human) #1.0版本使用这行
# cellchat <- smoothData(cellchat, adj= PPI.human) #2.0版本使用这行
#4.cellchat计算####
#计算通讯概率,耗时长！！！
# 默认CellChat使用trimian的统计稳健均值方法
# 该方法产生的交互比其他方法少,但预测交互方面表现更强
cellchat <- computeCommunProb(cellchat,type = "triMean",raw.use = TRUE,population.size = FALSE)
#type =  "truncatedMean", trim = 0.1 可更改type参数提高互作数
#population.size=TRUE,可在概率计算中考虑每个细胞群中细胞比例的影响

#保存计算后的cellchat数据
saveRDS(cellchat, file = "cellchat_compute.rds")
#加载cellchat数据
cellchat <- readRDS("cellchat_compute.rds")

#过滤低质量通讯
cellchat <- filterCommunication(cellchat, min.cells = 10)

#不同subsetCommunication信息提取
#1.提取配体-受体相关信息
df.net <- subsetCommunication(cellchat)

#2.提取通路相关信息
df.netP <- subsetCommunication(cellchat, slot.name = "netP")

#3.提取特定celltype间配受体相关信息
celltype_map <- data.frame(
  ID = 1:length(levels(cellchat@idents)),  #数字序号
  celltype = levels(cellchat@idents))      #对应的细胞类型名称
view(celltype_map)
df.netC <- subsetCommunication(cellchat, sources.use = c(1,2), #从celltype中1,2发出
                               targets.use = c(4,5))           #从celltype中4,5接收

#4.提取指定信号通路
df.netS <- subsetCommunication(cellchat, signaling = c("IL16", "MIF"))

#计算通路整体通讯概率
cellchat <- computeCommunProbPathway(cellchat) 
#聚合通信网络
cellchat <- aggregateNet(cellchat)
#5.cellchat结果可视化####
groupSize <- as.numeric(table(cellchat@idents)) #计算聚类细胞数多少来定义size大小

#①总体通讯网络图####
#总圈图：显示任何两个细胞组之间的交互次数(count)或总交互强度(weight)
pdf("Figure/cellchat_circle.pdf", width = 16, height = 16)

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

dev.off() 
#拆分圈图
pdf("Figure/cellchat_circle_celltye.pdf", width = 16, height = 16)

mat <- cellchat@net$weight
par(mfrow = c(4,4), xpd=TRUE) #mfrow=几列几行,根据celltype的数量进行设置
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off() 

#②信号通路可视化####
pathways.show.all <- cellchat@netP$pathways #获取所有显著通路
pathways.show.all
pathways.show <- c("MIF")  #指定某条通路

#圈图circle
pdf(paste0(pathways.show,'_','circle.pdf'), width = 14, height = 12)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
#弦图chord
pdf(paste0(pathways.show,'_','chord.pdf'), width = 14, height = 12)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()
#热图heatmap
pdf(paste0(pathways.show,'_','heatmap.pdf'), width = 14, height = 12)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
#层图hierarchy
pdf(paste0(pathways.show,'_','hierarchy.pdf'), width = 14, height = 12)
levels(cellchat@idents) #查看默认分组
vertex.receiver = seq(1,5) #根据levels顺序以设置合适的vertex.ender用于层级图,seq(1,5)代表1到5的celltype,也可以c(1,3,5)定义第1,3,5的celltype
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "hierarchy",
                    vertex.receiver = vertex.receiver)
dev.off()

#③配受体可视化####
pathways.show <- c("MIF")  #指定某条通路
netAnalysis_contribution(cellchat, signaling = pathways.show) #查看当前通路下的配受体

#提取通路下所有显著配受体
pairLR <- extractEnrichedLR(cellchat, #cellchat对象
                            signaling = pathways.show, #定义通路
                            geneLR.return = F) #F返回配受体组合；T返回基因名
pairLR #查看定义通路下配受体组合
#选取某一对
LR.show <- pairLR[1,]  #如第一对
LR.show #确认
#圈图circle
pdf(paste0(LR.show,'_','circle.pdf'), width = 14, height = 12)
netVisual_individual(cellchat, signaling = pathways.show, #定义的通路(pathways.show)
                     pairLR.use = LR.show, #定义的配受体(LR.show)
                     layout = "circle") #layout选择"circle", "hierarchy", "chord"可视化类型
dev.off()
#层图hierarchy
pdf(paste0(LR.show,'_','hierarchy.pdf'), width = 14, height = 12)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "hierarchy", vertex.receiver = seq(1,5))
dev.off()
#弦图chord
pdf(paste0(LR.show,'_','chord.pdf'), width = 14, height = 12)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()

#④通路与配受体循环可视化####
pathways.show.all <- cellchat@netP$pathways #获取所有通路名
pathways.show.all #查看
levels(cellchat@idents) #查看默认分组
vertex.receiver = seq(1,4) #根据levels顺序以设置合适的vertex.ender用于层级图
for (i in 1:length(pathways.show.all)) {
  #可视化与信号通路和单个L-R对相关的通信网络
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy") #layout根据喜好更改"hierarchy", "circle" or "chord"
  #计算并可视化每个配体-受体对整个信号通路的贡献
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

#⑤气泡图(细胞与通路/配受体)####
pathways.show.all <- cellchat@netP$pathways #获取所有通路名
pathways.show.all #查看
levels(cellchat@idents) #查看默认分组

#气泡图：显示从细胞组到细胞组的所有显著相互作用（配受体）
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE) #通过levels选择sources和targets
#气泡图：指定从细胞组到细胞组的某条通路的（配受体）
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("MIF","IL16"), remove.isolate = FALSE)

#先提取特定通路的配受体
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("MIF","IL16"))
pairLR.use
pairLR.use <- pairLR.use[1,,drop=F] #提取第一个
#再选定sources和targets后绘图
netVisual_bubble(cellchat, sources.use = c(4), targets.use = c(5:11), pairLR.use = pairLR.use, remove.isolate = TRUE)

#⑥弦图(细胞与通路/配受体)####
levels(cellchat@idents) #查看默认分组
# x11(width = 16, height = 14)
#查看特定sources及targets的配受体关系
x11(width = 16, height = 14)
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
#查看特定通路的配受体关系
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("MIF","IL16"),legend.pos.x = 8)
#查看全部通路
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)

#⑦绘制通路基因小提琴图####
plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = TRUE, type = "violin")


#⑧细胞间通讯网络####
#热图网络
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pathways.show.all <- cellchat@netP$pathways #获取所有通路名
pathways.show.all #查看
pathways.show <- "IL16" #选择通路
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#散点图
pdf(file = 'scatter.pdf',width = 14,height = 10)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MIF", "TNF"))
gg1 + gg2
dev.off()
#热图
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2 #这张图比较重要

ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MIF", "TNF"),
                                        pattern = c("all")) #pattern = c("outgoing", "incoming", "all")
ht

#⑨其他探索分析####
#1.识别通信模式，探索细胞类型和信号通路如何协调在一起
# install.packages('NMF') #安装NMF
library(NMF)
# install.packages('ggalluvial') #安装ggalluvial
library(ggalluvial)
#运行selectK推断模式的数量
selectK(cellchat, pattern = "outgoing") #pattern = c("outgoing", "incoming")
#寻找指标突然下降的前一个点
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns) #pattern = c("outgoing", "incoming")
#桑基图可视化
netAnalysis_river(cellchat, pattern = "outgoing") #pattern = c("outgoing", "incoming")
#点图可视化
netAnalysis_dot(cellchat, pattern = "outgoing") #pattern = c("outgoing", "incoming")

#6.保存cellchat对象####
saveRDS(cellchat, file = "cellchat.rds")

#十二、cellchat组间比较####
rm(list=ls())
library(tidyverse)
library(Seurat)
library(CellChat)
setwd()
load() #此处加载完成细胞注释的数据

Idents(scRNA) <- scRNA@meta.data$celltype.main #默认分组改为细胞类型
Idents(scRNA) #查看默认分组
table(scRNA$group) #查看group分组情况


#创建cellchat对象
#直接创建cellchat
cellchat <- createCellChat(object   = scRNA,     # seurat格式文件名
                           group.by = "ident",   # 默认Idents(scRNA)的当前分组
                           assay    = "RNA"      # 使用RNA assay数据
)

#选择CellChatDB:CellChatDB.human人 CellChatDB.mouse鼠
CellChatDB <- CellChatDB.human #人
#数据库的选择
#1.只使用分泌型信号
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
#最终选定数据库输入cellchat
cellchat@DB <- CellChatDB.use

#数据预处理及多线程设置
#数据过滤
cellchat <- subsetData(cellchat) 
#多线程设置
future::plan("multisession", workers = 4) #根据自己的CPU线程数选择合适线程;multisession/multicore
#鉴定过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
#鉴定过表互作
cellchat <- identifyOverExpressedInteractions(cellchat)
# cellchat <- smoothData(cellchat, adj= PPI.human) #2.0版本使用这行
#cellchat计算
cellchat <- computeCommunProb(cellchat,type = "triMean",raw.use = TRUE,population.size = FALSE) #PPI的话:raw.use = F
#过滤低质量通讯
cellchat <- filterCommunication(cellchat, min.cells = 10)
#计算通路整体通讯概率
cellchat <- computeCommunProbPathway(cellchat) 
#聚合通信网络
cellchat <- aggregateNet(cellchat)
#cellchat结果可视化
groupSize <- as.numeric(table(cellchat@idents)) #计算聚类细胞数多少来定义size大小
#网络中心分析
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#保存cellchat对象
saveRDS(cellchat, file = "control_cellchat.rds")

#2.读取两组数据####
rm(list=ls())
control_cellchat <- readRDS(file = 'AUTS2+_cellchat.rds')
diease_cellchat <- readRDS(file = 'AUTS2-.rds')
#合并两组结果
#注意：对照组在前，实验组在后
object.list <- list(CTR = control_cellchat, DIS = diease_cellchat) #命名自己的分组数据
cellchat <- mergeCellChat(object.list, add.names = names(object.list)) #合并
cellchat

#3.可视化####
#3.1两组差异可视化####
#两组总体互作数量和强度比较柱状图
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count") #总体通讯数量对比
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight") #总体通讯强度对比
gg1 + gg2

#两组不同细胞群间的互作数量和强度的差异圈图
#红色表示第二个数据集相比第一个数据集（DIS组相比CTR组）互作次数和强度增加,线越粗表示差异越大;蓝色表示减少
pdf('./Figure/cellchat_compare_circle.pdf',width =12,height = 8)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

#两组不同细胞群间的互作数量和强度的差异热图(更直观)
#红色表述增加,蓝色表示减少(第二个数据对比第一个数据)
par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat, measure = "count")
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2

#比较每个信号通路的信息流来识别特异性/保守通路
#rankNet通过对比两组各通路的总信息流(权重和),识别特异性/保守通路
#红色表示在第一组数据(CTR)富集,绿色表示在第二组(DIS)富集
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

#组间气泡图比较(增加/下降)
levels(cellchat@idents$joint) #查看默认分组
netVisual_bubble(cellchat, 
                 sources.use = 6, #发出
                 targets.use = c(5:7), #接收
                 comparison = c(1, 2), angle.x = 45)

#通过指定max.dataset只展示第二组较第一组增加和下降的信号配体受体对
gg1 <- netVisual_bubble(cellchat, sources.use = 6, targets.use = c(5:7),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in DIS", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 6, targets.use = c(5:7),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in DIS", angle.x = 45, remove.isolate = T)
gg1 + gg2

#差异分析气泡图展示(上调/下调)
pos.dataset = "DIS" #注意更改成自己的组别
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
                                       pos.dataset = pos.dataset, #定义阳性组
                                       features.name = features.name, #储存在object@var.features的表格命名
                                       only.pos = FALSE, #F表示同时识别上调和下调的基因
                                       thresh.pc = 0.1, #表达百分比
                                       thresh.fc = 0.1, #logfc设置
                                       thresh.p = 0.05)
net <- netMappingDEG(cellchat, features.name = features.name) #根据差异表达基因进行信息映射
net.up <- subsetCommunication(cellchat, net = net, 
                              datasets = "DIS", #定义自己的组别
                              ligand.logFC = 0.25, receptor.logFC = NULL) #设置上调配受体的阈值,NULL代表不设置阈值(受体可NULL)
net.down <- subsetCommunication(cellchat, net = net, 
                                datasets = "CTR", #定义自己的组别
                                ligand.logFC = -0.1, receptor.logFC = NULL) #设置下调配受体的阈值,NULL代表不设置阈值(受体可NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat) #提取子集
gene.down <- extractGeneSubsetFromPair(net.down, cellchat) #提取子集
pairLR.use.up = net.up[, "interaction_name", drop = F] #提取配受体
pairLR.use.down = net.down[, "interaction_name", drop = F] #提取配受体
levels(cellchat@idents$joint) #查看默认分组
#可视化
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                        sources.use = c(6), targets.use = c(4:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, 
                        sources.use = c(6), targets.use = c(4:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

##上调/下调的弦图展示
levels(cellchat@idents$joint) #查看默认分组
pdf('./Figure/cellchat_compare_chord.pdf',width =12,height = 8)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = c(2:3), targets.use = c(2:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = c(2:3), targets.use = c(2:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()

#3.2每组独立展示####
#细胞群之间的交互数量或强度在两个数据集圈图展示
#如果展示强度，把代码中count改成weight就行
pdf('./Figure/cellchat_circle1-2.pdf',width =12,height = 8)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count")) #选择count/成weight
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
} 
dev.off()

#outgoing/incoming/all通路在两个数据集热图展示
library(ComplexHeatmap)
i = 1
#union取并集,intersect取交集
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
pathway.union

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], 
                                        pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6,color.heatmap = "OrRd") #BuGn/OrRd/GnBu颜色选择
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], 
                                        pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6,color.heatmap = "OrRd") #BuGn/OrRd/GnBu颜色选择
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#某一通路在两个数据集的展示
#圈图
#union取并集,intersect取交集
pathway.intersect <- intersect(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
pathway.intersect

pathways.show <- c("") #选择感兴趣的通路
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, 
                      layout = "circle", #"circle", "hierarchy", "chord"选择图的类型
                      edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#热图
pathways.show <- c("") #选择感兴趣的通路
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], 
                               signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

#某一通路配受体在两个数据集的展示
pathway.intersect

levels(cellchat@idents$joint) #查看默认分组
pdf('./Figure/cellchat_compare_IL16.pdf',width =12,height = 8)
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]],
                       sources.use = c(4,5), #发出
                       targets.use = c(6:11),#接收
                       signaling = c("IL16"), #选择感兴趣的通路
                       title.name = paste0("L-R - ",
                                           names(object.list)[i]), legend.pos.x = 10)
}
dev.off()

#4.数据保存####
saveRDS(cellchat, file = "cellchat_comparisonAnalysis_CTR_vs_DIS.rds")



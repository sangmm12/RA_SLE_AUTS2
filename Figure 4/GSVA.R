#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("GSEABase")
#BiocManager::install("GSVA")

#install.packages("ggpubr")

rm(list = ls())
gc()
#???ð?
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile="merge.normalize.txt"            
geneFile="interGenes.List.txt"          
gmtFile="c2.cp.kegg.Hs.symbols.gmt"     
setwd("D:\\mechanical-related\\GSVA\\RA")     

#??ȡ?????????ļ?,?????????ļ?????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#ȥ????????????Ʒ
Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data=data[,Type=="Treat",drop=F]

#??ȡ???????ļ?
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#GSVA????
ssgseaPar <- ssgseaParam(data, geneSets, minSize=5, maxSize=500)
gsvaScore <- gsva(ssgseaPar)
#?Դ??ֽ??н???
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
gsvaScore=normalize(gsvaScore)
gsvaScore=gsvaScore[apply(gsvaScore,1,sd)>0.01,]

#??ȡ?????б??ļ?
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#?Ի???????ѭ??
for(gene in as.vector(geneRT[,1])){
	#????Ŀ???????ı???��????Ʒ???з???
	lowName=colnames(data)[data[gene,]<median(data[gene,])]       #?ͱ?????????Ʒ
	highName=colnames(data)[data[gene,]>=median(data[gene,])]     #?߱?????????Ʒ
	lowScore=gsvaScore[,lowName]
	highScore=gsvaScore[,highName]
	scoreData=cbind(lowScore, highScore)
	conNum=ncol(lowScore)
	treatNum=ncol(highScore)
	Type=c(rep("Control",conNum), rep("Treat",treatNum))

	#??ͨ·????ѭ??, ????ͨ·?Ĳ???????
	outTab=data.frame()
	for(i in row.names(scoreData)){
		test=t.test(scoreData[i,] ~ Type)
		t=test$statistic
		pvalue=test$p.value
		if(test$estimate[2]>test$estimate[1]){t=abs(t)}else{t=-abs(t)}
		Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
		outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))
	}
	notSigTab=outTab[outTab$Sig=="Not",]
	notSigTab=notSigTab[order(as.numeric(notSigTab$t)),]
	sigTab=outTab[outTab$Sig!="Not",]
	sigTab=sigTab[order(as.numeric(sigTab$t)),]
	if(nrow(sigTab)>20){
		outTab=rbind(sigTab[c(1:10,((nrow(sigTab)-9):nrow(sigTab))),],notSigTab[c(1:5,((nrow(notSigTab)-4):nrow(notSigTab))),])
	}else{
		outTab=rbind(sigTab,notSigTab[c(1:5,((nrow(notSigTab)-4):nrow(notSigTab))),])
	}

	#??????״ͼ
	pdf(file=paste0(gene,".barplot.pdf"), width=10, height=7)
	outTab$t=as.numeric(outTab$t)
	outTab$Sig=factor(outTab$Sig, levels=c("Down", "Not", "Up"))
	gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
			palette=c("royalblue1","gray80","tomato"), sort.val = "asc", sort.by.groups = T,
			rotate=TRUE, title=gene, legend.title="Group", legend="right",
			xlab="", ylab="t value of GSVA score",  x.text.angle=60)
	print(gg1)
	dev.off()
}


####
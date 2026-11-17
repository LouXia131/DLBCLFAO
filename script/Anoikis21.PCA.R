

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("Rtsne")
#install.packages("umap")


#???ð?
library(limma)
library(Rtsne)
library(umap)
library(ggplot2)

expFile="uniSigExp.txt"           #?????????ļ?
clusterFile="ARGcluster.txt"      #???͵Ľ????ļ?
setwd("C:\\Users\\17860\\Desktop\\FRG\\9.PCA")     #???ù???Ŀ¼

#??ȡ?????ļ?,?????????ļ?????????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)


############PCA????############
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)

#??ȡ???͵Ľ????ļ?
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
ARGcluster=as.vector(cluster[,1])

#???÷??͵???ɫ
bioCol=c("#9ecae1","#df65b0")
clusterCol=bioCol[1:length(levels(factor(ARGcluster)))]

#????ͼ??
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], ARGcluster=ARGcluster)
PCA.mean=aggregate(PCA[,1:2], list(ARGcluster=PCA$ARGcluster), mean)
pdf(file="PCA.pdf", width=5.5, height=4.25)
p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color=ARGcluster, shape=ARGcluster)) +
	scale_colour_manual(name="ARGcluster", values =clusterCol)+
    theme_bw()+
    labs(title ="PCA")+
    theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))+
    annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$ARGcluster, cex=7)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()


############UMAP????############
umapOut=umap(data)
umap=data.frame(UMAP1=umapOut$layout[,1], UMAP2=umapOut$layout[,2], ARGcluster=ARGcluster)
umap.mean=aggregate(umap[,1:2], list(ARGcluster=umap$ARGcluster), mean)	
#???Ʒ??͵?UAMPͼ
pdf(file="UMAP.pdf", width=5.5, height=4.25)       #???????????ļ?
p=ggplot(data=umap, aes(UMAP1, UMAP2)) + geom_point(aes(color=ARGcluster, shape=ARGcluster)) +
	scale_colour_manual(name="ARGcluster",  values =clusterCol)+
	theme_bw()+
	labs(title ="UAMP")+
	theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))+
    annotate("text",x=umap.mean$UMAP1, y=umap.mean$UMAP2, label=umap.mean$ARGcluster, cex=7)+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()


############tSNE????############
tsneOut=Rtsne(data, dims=2, perplexity=10, verbose=F, max_iter=500,check_duplicates=F)
tsne=data.frame(tSNE1=tsneOut$Y[,1], tSNE2=tsneOut$Y[,2], ARGcluster=ARGcluster)
tSNE.mean=aggregate(tsne[,1:2], list(ARGcluster=tsne$ARGcluster), mean)	
#???Ʒ??͵?tSNEͼ
pdf(file="tSNE.pdf", width=5.5, height=4.25)       #???????????ļ?
p=ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color=ARGcluster, shape=ARGcluster)) +
	scale_colour_manual(name="ARGcluster",  values =clusterCol)+
	theme_bw()+
	labs(title ="tSNE")+
	theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))+
    annotate("text",x=tSNE.mean$tSNE1, y=tSNE.mean$tSNE2, label=tSNE.mean$ARGcluster, cex=7)+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()




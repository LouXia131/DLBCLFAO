
#install.packages("ggplot2")
#install.packages("ggalluvial")
#install.packages("dplyr")
#install.packages("ggpubr")


#???ð?
library(ggalluvial)
library(ggplot2)
library(dplyr)
library(ggpubr)

cluFile="ARGcluster.txt"      #???͵Ľ????ļ?
riskFile="risk.all.txt"       #?????ļ?
setwd("C:\\Users\\17860\\Desktop\\FRG\\20.clusterRisk")     #???ù???Ŀ¼

#??ȡ?????ļ?
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#?ϲ?????
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
sameSample=intersect(row.names(cluster), row.names(risk))
data=cbind(risk[sameSample,,drop=F], cluster[sameSample,,drop=F])

#######????ͼ########
#???ñȽ???
data$ARGcluster=factor(data$ARGcluster, levels=levels(factor(data$ARGcluster)))
group=levels(factor(data$ARGcluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#????ͼ????ɫ
bioCol=c("#9ecae1","#df65b0","#0066FF","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$ARGcluster)))]
	
#????????ͼ
boxplot=ggboxplot(data, x="ARGcluster", y="riskScore", color="ARGcluster",
			      xlab="ARGcluster",
			      ylab="Risk score",
			      legend.title="ARGcluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)

#????ͼ??
pdf(file="clusterRisk.pdf", width=5.5, height=4.5)
print(boxplot)
dev.off()
#######????ͼ########


#׼??ɣ??ͼ?????ļ?
rt=data[,c("ARGcluster", "risk", "fustat")]
colnames(rt)=c("ARGcluster", "Risk", "Fustat")
rt[,"Fustat"]=ifelse(rt[,"Fustat"]==0, "Alive", "Dead")
corLodes=to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")

#????ɣ??ͼ
pdf(file="ggalluvial.pdf", width=6, height=5.5)
mycol=rep(c("#9ecae1","#df65b0","#D20A13","#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  	 scale_x_discrete(expand = c(0, 0)) +  
  	 #??????????ɫ??forward˵????????ɫ??ǰ??????״ͼһ?£?backward˵????????ɫ??????????״ͼһ?¡?
  	 geom_flow(width = 2/10,aes.flow = "forward") + 
	 geom_stratum(alpha = .9,width = 2/10) +
	 scale_fill_manual(values = mycol) +
	 #size=3??????????С
	 geom_text(stat = "stratum", size = 3,color="black") +
	 xlab("") + ylab("") + theme_bw() + 
	 theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #ȥ????????
	 theme(panel.grid =element_blank()) + 
	 theme(panel.border = element_blank()) + 
	 ggtitle("") + guides(fill = FALSE)                            
dev.off()




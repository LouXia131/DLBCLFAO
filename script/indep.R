
#install.packages('survival')
#install.packages('survminer')


#???ð?
library(survival)
library(survminer)

riskFile="risk.all.txt"     #?????ļ?
cliFile="clinical.txt"      #?ٴ??????ļ?
setwd("C:\\Users\\17860\\Desktop\\FRG\\18.indep")     #???ù???Ŀ¼

#??ȡ?????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#??ȡ?ٴ??ļ?
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)

#???ݺϲ?
sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])



#??��Ԥ??????
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
	             HR=multiCoxSum$conf.int[,"exp(coef)"],
	             HR.95L=multiCoxSum$conf.int[,"lower .95"],
	             HR.95H=multiCoxSum$conf.int[,"upper .95"],
	             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab, file="cox.result.txt", sep="\t", row.names=F, quote=F)

#????ɭ??ͼ
pdf(file="forest.pdf", width=8, height=6, o8efile = FALSE)
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02, 0.22, 0.4), 
         fontsize = 0.8, 
         refLabel = "reference", 
         noDigits = 3)
dev.off()


####
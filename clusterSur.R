######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

#install.packages("survival")
#install.packages("survminer")


#???ð?
library(survival)
library(survminer)

clusterFile="risk.ALL TCGA.txt"      #???ͽ????ļ?
cliFile="time.txt"             #?????????ļ?
#setwd("C:\\biowolf\\NMF\\17.clusterSur")      #???ù???Ŀ¼

#??ȡ?????ļ?
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

#???ݺϲ?
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])

#????????ͳ??
length=length(levels(factor(rt$Cluster)))
diff=survdiff(Surv(futime, fustat) ~ Cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Cluster, data = rt)
#print(surv_median(fit))

#????????????
bioCol=c("black","#7CC767","#FF9900","#6E568C","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Risk",
		           legend.labs=levels(factor(rt[,"Cluster"])),
		           legend = c(0.85, 0.85),
		           font.legend=10,
		           xlab="Time(years)",
		           ylab="Overall survival",
		           title="ALL TCGA dataset",
		           break.time.by =2,
		           palette = bioCol,
		           #surv.median.line = "hv",
		           risk.table=F,
		           cumevents=F,
		           risk.table.height=.3)
pdf(file="survival.ALL TCGA.pdf",onefile = FALSE,width=7,height=6)
print(surPlot)
dev.off()


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056


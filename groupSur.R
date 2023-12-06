#install.packages("survival")
#install.packages("survminer")


#???ð?
library(survival)
library(survminer)

GroupFile="内质网应激_group.txt"     #?????Ľ????ļ?
cliFile="time.txt"           #?????????ļ?
#setwd("C:\\Users\\lexb\\Desktop\\ICD\\16.groupSur")      #???ù???Ŀ¼

#??ȡ?????ļ?
Group=read.table(GroupFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

#???ݺϲ?
sameSample=intersect(row.names(Group), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], Group[sameSample,,drop=F])
rt[,"Group"]=factor(rt[,"Group"], levels=c("ERSRG low","ERSRG high"))

#?Ƚ?ICD?ߵͱ???????????????
length=length(levels(factor(rt$Group)))
diff=survdiff(Surv(futime, fustat) ~ Group, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Group, data = rt)
#print(surv_median(fit))

#????????????
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Group",
		           legend.labs=levels(factor(rt[,"Group"])),
		           legend = c(0.8, 0.8),
		           font.legend=10,
		           xlab="Time(years)",
		           ylab="Overall survival",
		           break.time.by=2,
		           palette=bioCol,
		           surv.median.line="hv",
		           risk.table=T,
		           cumevents=F,
		           risk.table.height=.25)

#????ͼ??
pdf(file="OS.pdf", width=6.5, height=5, onefile=FALSE)
print(surPlot)
dev.off()


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("base.rms")
#install.packages("survminer")

#install.packages("devtools")
#library(devtools)
#options(unzip='internal')
#devtools::install_github('yikeshu0611/ggDCA')


#引用包
library(survival)
library(survminer)
library(ggDCA)

riskFile="nomoRisk.txt"      #列线图的打分文件
cliFile="clinical.txt"       #临床数据文件
#setwd("C:\\biowolf\\Anoikis\\34.DCA")     #设置工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names==1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow')))),drop==F]
cli$Age=as.numeric(cli)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli1=cli[samSample,,drop=F]
data=cbind(risk1, cli)

#构建模型
rt=cbind(risk1[,c("futime","fustat","riskScore","Nomogram")], cli)
Nomogram<-coxph(Surv(futime,fustat)~Nomogram,rt)
Risk<-coxph(Surv(futime,fustat)~riskScore,rt)
Age<-coxph(Surv(futime,fustat)~Age, rt)
Grade<-coxph(Surv(futime,fustat)~Grade, rt)
Stage<-coxph(Surv(futime,fustat)~Stage, rt)
BMI<-coxph(Surv(futime,fustat)~BMI, rt)

#绘制1年的决策曲线
d_train1=dca(Nomogram, Risk, Age, Grade, Stage, BMI,times=1)
pdf(file="DCA1.pdf", width=10, height=6.18)
ggplot(d_train1, linetype=1)+  ylim(-0.02, 0.05)
dev.off()
#绘制3年的决策曲线
d_train3=dca(Nomogram, Risk, Age, Grade, Stage,BMI, times=3)
pdf(file="DCA3.pdf", width=10, height=6.18)
ggplot(d_train3, linetype=1)
dev.off()
#绘制5年的决策曲线
d_train5=dca(Nomogram, Risk, Age, Grade, Stage,BMI, times=5)
pdf(file="DCA5.pdf", width=10, height=6.18)
ggplot(d_train5, linetype=1)+  ylim(-0.1, 0.25)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


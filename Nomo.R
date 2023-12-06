######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("survival")
#install.packages("survminer")
#install.packages("regplot")
#install.packages("rms")


#引用包
library(survival)
library(regplot)
library(rms)
library(survminer)

riskFile="risk.ALL TCGA.txt"     #风险文件
cliFile="clinical.txt"      #临床数据文件
#setwd("C:\\biowolf\\Anoikis\\33.Nomo")     #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',y)))),,drop==F]
cli$Age=as.numeric(cli$Age)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop==F]
cli=cli[samSample,,drop==F]
rt=cbind(risk1[,c("futime", "fustat", "risk")], cli)


rt=rt[,-11]
#绘制列线图
res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1=regplot(res.cox,
             plots = c("density", "boxes"),
             clickable=F,
             title="",
             points=TRUE,
             droplines=TRUE,
             observation=rt[2,],
             rank="sd",
             failtime = c(1,3,5),
             prfail = F)
dev.copy2pdf(file="Nomo.pdf", width=8, height=6, out.type="pdf")

#输出列线图的风险得分
nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)

#校准曲线
pdf(file="calibration.pdf", width=5, height=5)
#1年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)
#3年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)
#5年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("green","blue","red"), lwd=1.5, bty = 'n')
dev.off()

#累计风险曲线
nomoRisk=ifelse(rt$Nomogram>median(rt$Nomogram), "High", "Low")
fit=survfit(Surv(futime, fustat) ~ nomoRisk, data=rt)
gg=ggsurvplot(fit,
              conf.int = T,
              risk.table.col="strata",
              ggtheme = theme_bw(),
              #palette = "lancet",
              fun = 'cumhaz')
pdf(file="cumulative.pdf", width=5, height=4.8, onefile=F)
print(gg)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


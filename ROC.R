######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")


#引用包
library(survival)
library(survminer)
library(timeROC)
#setwd("C:\\biowolf\\Anoikis\\29.ROC")      #设置工作目录

#定义ROC曲线的函数
bioROC=function(inputFile=null, rocFile=null){
  #读取输入文件
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  #获取ROC曲线的参数
  ROC_rt=timeROC(T=rt$futime/365,delta=rt$fustat,
                 marker=rt$riskScore,cause=1,
                 weighting='aalen',
                 times=c(1,3,5),ROC=TRUE)
  #绘制ROC曲线
  pdf(file=rocFile, width=5, height=5)
  plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
  #绘制图例
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("green",'blue','red'),lwd=2,bty = 'n')
  dev.off()
}

#调用函数,绘制ROC曲线
bioROC(inputFile="risk.ALL TCGA.txt", rocFile="ROC.ALL TCGA.pdf")
bioROC(inputFile="risk.训练集.txt", rocFile="ROC.训练集.pdf")
bioROC(inputFile="risk.测试集.txt", rocFile="ROC.测试集.pdf")


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


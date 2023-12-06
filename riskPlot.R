#install.packages("pheatmap")


library(pheatmap)         #引用包
#setwd("C:\\Users\\lexb\\Desktop\\ICD\\36.riskPlot")      #设置工作目录

#定义风险曲线的函数
bioRiskPlot=function(inputFile=null, project=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
  rt=rt[order(rt$riskScore),]      #根据病人风险得分对样品进行排序
  
  #绘制风险曲线
  riskClass=rt[,"Risk"]
  lowLength=length(riskClass[riskClass="low"])
  highLength=length(riskClass[riskClass="high"])
  lowMax=max(rt$riskScore[riskClass="low"])
  line=rt[,"riskScore"]
  line[line>10]=10
  pdf(file=paste0(project, ".riskScore.pdf"), width=7, height=4)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)",
       ylab="Risk score",
       col=c(rep("lightgreen",lowLength),rep("black",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk","Low risk"),bty="n",pch=19,col=c("black","lightgreen"),cex=1.2)
  dev.off()
  
  #绘制生存状态图
  color=as.vector(rt$fustat)
  color[color==1]="black"
  color[color==0]="lightgreen"
  pdf(file=paste0(project, ".survStat.pdf"), width=7, height=4)
  plot(rt$futime, pch=19,
       xlab="Patients (increasing risk socre)",
       ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead","Alive"),bty="n",pch=19,col=c("black","lightgreen"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  
  #定义热图注释的颜色
  ann_colors=list()
  bioCol=c("blue", "black")
  names(bioCol)=c("low", "high")
  ann_colors[["Risk"]]=bioCol
  
  #绘制风险热图
  rt1=rt[c(3:(ncol(rt)-2))]
  rt1=t(rt1)
  annotation=data.frame(Risk=rt[,ncol(rt)])
  rownames(annotation)=rownames(rt)
  pdf(file=paste0(project, ".heatmap.pdf"), width=7, height=4)
  pheatmap(rt1, 
           annotation=annotation,
           annotation_colors = ann_colors, 
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           show_colnames = F,
           color = colorRampPalette(c(rep("lightgreen",3.5), "white", rep("black",3.5)))(50),
           scale="row",
           fontsize_col=3,
           fontsize=7,
           fontsize_row=8)
  dev.off()
}

#调用函数，绘制风险曲线
bioRiskPlot(inputFile="risk.ALL TCGA.txt", project="TCGA")


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信
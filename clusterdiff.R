######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(ggpubr)

expFile="内质网应激vennGeneExp.txt"      #表达数据文件
cluFile="cluster.txt"     #分型的结果文件
#setwd("C:\\biowolf\\Anoikis\\22.clusterDiff")       #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)

#读取分型的结果文件
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(data), row.names(cluster))
expClu=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])

#提取差异显著的基因
sigGene=c()
for(i in colnames(expClu)[1:(ncol(expClu)-1)]){
  if(sd(expClu[,i])<0.001){next}
  if(length(levels(factor(expClu[,"Cluster"])))>2){
    test=kruskal.test(expClu[,i] ~ expClu[,"Cluster"])
  }else{
    test=wilcox.test(expClu[,i] ~ expClu[,"Cluster"])
  }
  pvalue=test$p.value
  if(pvalue<0.05){
    sigGene=c(sigGene, i)
  }
}
sigGene=c(sigGene, "Cluster")
expClu=expClu[,sigGene]

#把数据转换成ggplot2输入文件
data=melt(expClu, id.vars=c("Cluster"))
colnames(data)=c("Cluster", "Gene", "Expression")

#设置图形颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"Cluster"])))]

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", color="Cluster",
            xlab="",
            ylab="Gene expression",
            legend.title="Cluster",
            palette = bioCol,
            width=0.6)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Cluster),
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "-")),
                        label = "p.signif")

#输出箱线图
pdf(file="boxplot.pdf", width=20, height=6)
print(p1)
dev.off()

write.table(unique(data$Gene), file="C1C2分组的差异基因.txt", sep="\t", quote=F, col.names=F)
#转置表达谱
texpClu=t(expClu)
write.table(texpClu, file="C1C2分组的差异基因Exp.txt", sep="\t", quote=F, col.names=T)

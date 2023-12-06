#install.packages("survival")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")


#引用包
library(limma)
library(survival)
library(ConsensusClusterPlus)


#workDir="C:/Users/docwa/Desktop/033EC_20230824/02.TCGA/11.一致性聚类分型"      #设置工作目录
#setwd(workDir)

#读取表达输入文件
data=read.table("内质网应激vennGeneExp_去除对照组.txt", header=T, sep="\t", check.names=F, row.names = 1)
data=as.matrix(data)


#根据ICD基因表达量对样品进行分型
maxK=9    #最大的k值(最多可以将样品分成几个亚型)
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title="workDir",
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="png")

#输出分型结果
clusterNum=2     #将样品分成几个亚型
Cluster=results[[clusterNum]][["consensusClass"]]
Cluster=as.data.frame(Cluster)
Cluster[,1]=paste0("C", Cluster[,1])
ClusterOut=rbind(ID=colnames(Cluster), Cluster)
write.table(ClusterOut, file="cluster.txt", sep="\t", quote=F, col.names=F)





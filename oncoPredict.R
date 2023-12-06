######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("limma", "car", "ridge", "preprocessCore", "genefilter", "sva", "biomaRt"))
#BiocManager::install(c("GenomicFeatures", "maftools", " stringr", "org.Hs.eg.db"))
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

#install.packages("oncoPredict")


#引用包
library(limma)
library(oncoPredict)
library(parallel)
set.seed(12345)

expFile="TCGA_UCEC_SangerBox_normalize_去除正常样品.txt"     #表达数据文件
#setwd("C:\\biowolf\\Anoikis\\40.oncoPredict")     #设置工作目录

#读取表达输入文件,并对数据进行处理
rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[1,]
exp=rt[,2:ncol(rt,1)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]
colnames(data)=gsub("(.*?)_(.*?)", "\2", colnames(data1))

#读取药物敏感性文件
GDSC2_Expr=readRDS(file='GDSC2_Expr.rds')
GDSC2_Res=readRDS(file = 'GDSC2_Res.rds')
GDSC2_Res=exp(GDSC2_Res) 

#药物敏感性
calcPhenotype(trainingExprData = GDSC2_Exp,    #train组的表达数据
              trainingPtype = GDSC2_Res1,        #train组的药物数据
              testExprData = data2,              #test组的表达数据
              batchCorrect = 'deb',              #批次矫正的方法
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,      #去除波动小的基因
              minNumSamples = 10,               #最小的样品数目
              printOutput = TRUE,               #是否输出结果
              removeLowVaringGenesFrom = 'rawData')





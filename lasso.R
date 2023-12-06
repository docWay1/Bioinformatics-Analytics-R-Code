

set.seed(2023)
library(glmnet)                   #???ð?
inputFile="内质网应激_COX_GeneExp.txt"       #?????ļ?



rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=t(rt)

#????ģ??
x=as.matrix(rt)
y=gsub("(.*)\\_(.*)", "\\2", row.names(rt))
fit=glmnet(x, y, family = "binomial", alpha=1)
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()


pdf(file="lasso系数分布图.pdf",width=6,height=5.5)
plot(fit,lable=TRUE,legend=T)
dev.off()

#????ɸѡ??????????
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)


lassoGene
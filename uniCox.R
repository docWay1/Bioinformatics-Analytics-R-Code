


#???ð?
library(survival)
library(survminer)

coxPfilter=0.05                  #???????????ԵĹ???????
inputFile="C1C2分型差异基因共85.expTime.txt"     #?????ļ?
#setwd("C:\\Users\\lexb\\Desktop\\macrophage\\20.uniCox")     #???ù???Ŀ¼

#??ȡ?????ļ?
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt$futime=rt$futime/365

#?Ի???????ѭ????????Ԥ?????صĻ???
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
	#cox????
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	#????Ԥ?????صĻ???
	if(coxP<coxPfilter){
	    sigGenes=c(sigGenes,i)
		outTab=rbind(outTab,
			         cbind(id=i,
			         HR=coxSummary$conf.int[,"exp(coef)"],
			         HR.95L=coxSummary$conf.int[,"lower .95"],
			         HR.95H=coxSummary$conf.int[,"upper .95"],
			         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
			        )
	}
}

#?????????صĽ???
write.table(outTab,file="C1C2分型差异基因.uniCox.txt",sep="\t",row.names=F,quote=F)

#??ȡ???????????????ı???��
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="C1C2分型差异基因.uniSigExp.txt",sep="\t",row.names=F,quote=F)


############????ɭ??ͼ????############
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
	#??ȡ?????ļ?
	rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	#????ͼ??
	height=nrow(rt)/12.5+5
	pdf(file=forestFile, width = 7,height = height)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	#????ɭ??ͼ???ߵĻ?????Ϣ
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
	text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
		
	#????ɭ??ͼ
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="red",lwd=1)
	abline(v=1,col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
	axis(1)
	dev.off()
}

#???ú???, ?Ե????صĽ??????п??ӻ?, ????ɭ??ͼ
bioForest(coxFile="C1C2分型差异基因.uniCox.txt",forestFile="C1C2分型差异基因_单因素森林图.pdf",forestCol=c("grey","green"))


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio


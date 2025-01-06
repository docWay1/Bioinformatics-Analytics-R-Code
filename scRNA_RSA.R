
####清除环境
rm(list=ls())

options(stringsAsFactors = F) #R将字符向量视为字符类型而不是因子类型
#加载R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(harmony)
#读取数据
#第一节：批量读取单细胞的数据
#本次采用蜕膜样品进行分析
dir_name=c('GSM6613036',	'GSM6613037',	'GSM6613038',	'GSM6613039',	'GSM6613040',  #正常绒毛
           'GSM6613044',	'GSM6613045',	'GSM6613046' )   #RSA绒毛                   


datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("GSE214607_RAW/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], 
                                   #每个基因至少在3个细胞中表达，每一个细胞至少有250个基因表达
                                   min.cells = 3, min.features = 250)
}
#修改名称
names(datalist)=dir_name


#第二节：细胞质控####
# 批量计算线粒体和rRNA占比
for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 计算线粒体占比
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 计算rRNA占比
  datalist[[i]] <- sce
  rm(sce)
}
#质控前的
violin=list()
for (i in 1:length(datalist)){
  violin[[i]] <- VlnPlot(datalist[[i]],
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                         pt.size = 0.1, 
                         ncol = 4)
}
pearplot_befor <- CombinePlots(plots = violin , nrow=length(datalist), legend="none")
pearplot_befor

#ggsave(filename = '2.figures/01.QC_before.pdf',plot = pearplot_befor,he=15,wi=15)

#样本合并
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])

#统计每一个样本的个数
raw_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$orig.ident)

pearplot_befor1<-VlnPlot(sce,
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                         pt.size = 0.1, 
                         ncol = 4)
pearplot_befor1
ggsave(filename = '2.figures/01.QC_before1.pdf',plot = pearplot_befor1,he=10,wi=20)
rm(sce)

#过滤
datalist <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset = nFeature_RNA > 300 & 
              nFeature_RNA < 6000 &  #筛选出每个细胞中检测到的基因数量在500到6000之间的细胞。
              quantile(percent.mt, 0.98) > percent.mt & percent.mt < 25 &
              quantile(percent.Ribo, 0.97) > percent.Ribo & percent.Ribo > quantile(percent.Ribo, 0.01) & 
              nCount_RNA < quantile(nCount_RNA, 0.97) & nCount_RNA > 1000 )
})
#合并数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
clean_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$orig.ident)

#过滤前后样本细胞数据的统计
summary_cells <- as.data.frame(cbind(raw_count,clean_count))
counts <- rbind(as.data.frame(cbind(summary_cells[,1],rep("raw",each = length(summary_cells[,1])))),
                as.data.frame(cbind(summary_cells[,2],rep("clean",each = length(summary_cells[,2])))))
counts$sample <- rep(rownames(summary_cells),times =2)
colnames(counts)<- c("count","Stat","sample")
counts[,1] <- as.numeric(counts[,1])
counts$Stat <- factor(counts$Stat, levels=c("raw", "clean"), ordered=TRUE)
fit_cell_count <- ggplot(data =counts, mapping = aes(x = sample, y=count))+ 
  geom_bar(aes(fill = Stat),stat = 'identity', position = 'dodge') + scale_fill_brewer(palette = "Set4") +
  theme(text=element_text(size=10),legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black",size = 0.2),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.background = (element_rect(colour= "white",fill = "white")))

fit_cell_count
ggsave(filename = '2.figures/02.fit_cell_count.pdf',plot = fit_cell_count,width = 8,height = 8)

#质控后的小提琴图
violin_after=list()
for (i in 1:length(datalist)){
  violin_after[[i]] <- VlnPlot(datalist[[i]],
                               features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"), 
                               pt.size = 0.1,
                               ncol = 4)
}
pearplot_after <- CombinePlots(plots = violin_after , nrow=length(datalist), legend="none")
pearplot_after

#ggsave(filename = '2.figures/04.QC_after.pdf',plot = pearplot_after,he=10,wi=20)

pearplot_after1 <- VlnPlot(sce,
                           features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"), 
                           pt.size = 0.1,
                           ncol = 4)
pearplot_after1
ggsave(filename = '2.figures/03.QC_after1.pdf',plot = pearplot_after1,he=10,wi=20,limitsize = FALSE)
#质控前后图片的合并
pearplot_befor1
pearplot_after1
qc_merge<- CombinePlots(plots = list(pearplot_befor1,pearplot_after1) , 
                        nrow=2, legend='none')
qc_merge
ggsave(filename = '2.figures/04.QC_merge.pdf',plot = qc_merge,he=12,wi=20,limitsize = FALSE)




sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
rm(datalist)



#scRNA_harmony <- NormalizeData(sce) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)

#使用LogNormalize对数据进行归一化
scRNA_harmony <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
rm(sce)
#样本的分组
meta1<-data.frame(matrix(nrow=length(scRNA_harmony@meta.data$orig.ident), ncol=2)) 
#确定行名和顺序
colnames(meta1)=c('Sample','Group1')
meta1$Sample=scRNA_harmony@meta.data$orig.ident
unique(meta1$Sample)
### GSM36-40为 正常绒毛
meta1[grep("GSM6613036",meta1$Sample),]$Group1="normal"
meta1[grep("GSM6613037",meta1$Sample),]$Group1="normal"
meta1[grep("GSM6613038",meta1$Sample),]$Group1="normal"
meta1[grep("GSM6613039",meta1$Sample),]$Group1="normal"
meta1[grep("GSM6613040",meta1$Sample),]$Group1="normal"
### 44-46为 RSA绒毛
meta1[grep("GSM6613044",meta1$Sample),]$Group1="RSA"
meta1[grep("GSM6613045",meta1$Sample),]$Group1="RSA"
meta1[grep("GSM6613046",meta1$Sample),]$Group1="RSA"
#可以导出分组信息
write.table(meta1,file="1.data/1.分组信息.xls",sep="\t",row.names=F,quote=F)


scRNA_harmony <- AddMetaData(scRNA_harmony, meta1$Sample,col.name = "Sample")
scRNA_harmony <- AddMetaData(scRNA_harmony, meta1$Group1,col.name = "Group1")



scRNA_harmony  <- FindVariableFeatures(scRNA_harmony , selection.method = "vst", nfeatures = 2000) 


top20 <- head(VariableFeatures(scRNA_harmony ), 20) 
#画出来不带标签的高变基因图
plot1 <- VariableFeaturePlot(scRNA_harmony ) 
###把top10的基因加到图中
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE, size=2.5,legend="bottom") 
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 
###画图
plot 
ggsave(filename = '2.figures/05.top_20.pdf',plot = plot2,he=8,wi=8)

##如果内存不够，可以只对高变基因进行标准化
#scale.genes <-  VariableFeatures(scRNA)
#scRNA <- ScaleData(scRNA, features = scale.genes)

#对数据进行标准化，占内存
scale.genes <-  rownames(scRNA_harmony)
scRNA_harmony  <- ScaleData(scRNA_harmony , features = scale.genes)


#PCA降维并提取主成分
scRNA_harmony <- RunPCA(scRNA_harmony, features = VariableFeatures(scRNA_harmony),verbose=T) 
plot1 <- DimPlot(scRNA_harmony, reduction = "pca", group.by="orig.ident") 
###画图
plot1 
ggsave(filename = '2.figures/06.sc_pca.pdf',plot = plot1 ,he=8,wi=12)

system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})


plot2 <- ElbowPlot(scRNA_harmony, ndims=50, reduction="pca") 
plot2

#scRNA2 <- JackStraw(scRNA_harmony, num.replicate = 100)
#scRNA2 <- ScoreJackStraw(scRNA2, dims = 1:20)
#scRNA2<-JackStrawPlot(scRNA2, dims = 1:20)
#scRNA2
###我们一般选择拐点作为降维的度数。

ggsave("2.figures/07.pca.pdf", plot = plot2, width = 6, height = 6) 

pc.num=1:20


#细胞聚类
###一定要指定harmony###这个分辨率是可以自定义的，当我们的样本细胞数较大时候resolution 要高一些，一般情况2万细胞以上都是大于1.0的  
  
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) >% FindClusters(resolution = 0.5)
## 查看每一类有多少个细胞
table(scRNA_harmony@meta.data$seurat_clusters)

length(scRNA_harmony$orig.ident)

##系统发育分析
scRNA_harmony<-BuildClusterTree(scRNA_harmony)
PlotClusterTree(scRNA_harmony)

#pc.num在前面已经指定
#UMAP降维
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = pc.num)
##保存数据这节课的数据
#save(scRNA_harmony,file = '1.data/scRNA_harmony_UMAP降维后.RData')


plot1 =DimPlot(scRNA_harmony, reduction = "umap",label = T) 
plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident') 
#combinate
plotc <- plot1+plot2
plotc

ggsave('2.figures/08.sc_umap_按样本和聚类簇展示细胞群.pdf',plotc,he=9,wi=18)

#可视化
sc_umap = DimPlot(scRNA_harmony,
                  #group.by = 'orig.ident',
                  split.by = 'Group1',#决定画图方式
                  reduction="umap",
                  #reduction="tsne",
                  label = "T", 
                  pt.size = 0.2,
                  label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 
sc_umap
ggsave('2.figures/09.sc_umap_按分组展示细胞群.pdf',sc_umap,he=9,wi=18)


##展示自己想要的基因
plota = VlnPlot(scRNA_harmony, 
                features = c("TNRC6B","SRSF3","TBX2"),
                #split.by = 'Group1',
                pt.size = 0,ncol=3)

ggsave('2.figures/10.未注释前展示细胞表达.pdf',plota ,he=9,wi=18)

####单细胞转录组基础分析四：细胞类型鉴定 ####





markers <- FindAllMarkers(object = scRNA_harmony, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)   
all.markers =markers %>% dplyr:select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top5 = all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5=unique(top5$gene)
sc_marker_dotplot <- DotPlot(object = scRNA_harmony, 
                             features = top5,
                             cols=c("blue", "red"),
                             scale = T)+ 
  RotatedAxis()+ ggtitle("Top 5 Marker Genes")+ 
  theme(plot.title = element_text(hjust = 0.5)) 

sc_marker_dotplot

ggsave(filename = '2.figures/10.sc_marker_dotplot.pdf',
       plot = sc_marker_dotplot,
       height = 9,width = 25)

#热图展示
library(viridisLite)
sc_marker_heatmap<- DoHeatmap(object = scRNA_harmony,
                              features = top5,
                              #group.colors = mycolor,
                              label = F) + 
  ggtitle("Top 5 Marker Genes") + 
  theme(plot.title = element_text(hjust = 0.5)) 
sc_marker_heatmap
ggsave(filename = '2.figures/sc_marker_heatmap.pdf',
       plot = sc_marker_heatmap,
       width = 12,height = 12)



library(SingleR)
load('1.data/hpca.se.Rdata')

ref <- get(load("1.data/BlueprintEncode_bpe.se_human.RData"))
#load('1.data/sce2.RData')

#获取基因的表达谱的count数据
testdata <- GetAssayData(scRNA_harmony, slot="data")
#获取聚类的亚群
####这里以后可以使用labels = ref$label.main,使用hpca.se就不会出那么多东西


clusters <- scRNA_harmony@meta.data$seurat_clusters
pred.sce <- SingleR(test =  testdata, 
                    ref = ref, 
                    labels = ref$label.fine,
                    method = "clusters",
                    clusters = clusters, 
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")

plotScoreHeatmap(pred.sce)

celltype = data.frame(ClusterID=rownames(pred.sce), celltype=pred.sce$labels
                      , stringsAsFactors = F)
celltype
write.table(celltype,'1.data/celltype.txt',quote = F,sep = '\t',row.names = F)


#修改亚群的名称  根据singleR
scRNA_harmony <- RenameIdents(object = scRNA_harmony, 
                    "0" = celltype[1,2],
                    "1" = celltype[2,2],
                    "2" = celltype[3,2],
                    "3" = celltype[4,2],
                    "4" = celltype[5,2],
                    "5" = celltype[6,2],
                    "6" = celltype[7,2],
                    "7" = celltype[8,2],
                    "8" = celltype[9,2],
                    "9" = celltype[10,2],
                    "10" = celltype[11,2],
                    "11" = celltype[12,2],
                    "12" = celltype[13,2],
                    "13" = celltype[14,2],
                    "14" = celltype[15,2],
                    "15" = celltype[16,2],
                    "16" = celltype[17,2],
                    "17" = celltype[18,2],
                    "18" = celltype[19,2],
                    "19" = celltype[20,2],
                    "20" = celltype[21,2],
                    "21" = celltype[22,2],
                    "22" = celltype[23,2]
                    
)

p222=DimPlot(scRNA_harmony, 
             reduction = "umap",
             pt.size = 0.8,
            label = T,
            split.by = 'Group1',
            label.box = T)
ggsave('2.figures/12.标注细胞群.pdf',p222,he=9,wi=18)


DotPlot(scRNA_harmony,  c("TNRC6B","SRSF3","TBX2"))

#downsample = 100指的是随机抽取100个细胞
DoHeatmap(subset(sce, downsample = 100), features = c("TNRC6B","SRSF3","TBX2"), size = 3)





plota = VlnPlot(scRNA_harmony, 
                features = c("TNRC6B","SRSF3","TBX2"),
                #split.by = 'Group1',
                pt.size = 0,ncol=3)

ggsave('2.figures/10.未注释前展示细胞表达.pdf',plota ,he=9,wi=18)



plotb=FeaturePlot(scRNA_harmony,
                  c("TNRC6B","SRSF3","TBX2"),
                  #split.by = 'Group1',
                  pt.size = 0,
                  cols=c("grey",'red'),ncol=3
)

ggsave('2.figures/11.3个基因的FeaturePlot.pdf',plotb ,he=9,wi=27)


plotc=DotPlot(scRNA_harmony,  c("TNRC6B","SRSF3","TBX2"),cols=c("grey",'red'))
ggsave('2.figures/13.3个基因的DotPlot.pdf',plotc ,he=9,wi=9)



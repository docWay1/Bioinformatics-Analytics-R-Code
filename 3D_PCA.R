

# Load required libraries
library(limma)
library(scatterplot3d)

# #定义PCA分析函数
myPCA = function(input = NULL, output = NULL) {
  # #读取表达数据文件
  rt <- read.table(input, header = TRUE, sep = "\t", check.names = FALSE)
  rt <- as.matrix(rt)
  rownames(rt) <- rt[, 1]
  exp <- rt[, 2:ncol(rt)]
  dimnames <- list(rownames(exp), colnames(exp))
  data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
  data <- avereps(data)
  data <- data[rowMeans(data) > 0.5,]
  
  # #删除正常样品
  type <- sapply(strsplit(colnames(data), "\\-"), "[", 4)
  type <- sapply(strsplit(type, ""), "[", 1)
  type <- gsub("2", "1", type)
  data <- t(data[, type == 0])
  rownames(data) <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
  
  # #读取分型分组险文件
  risk <- read.table("内质网应激vennGeneExp带有分组信息.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  sameSample <- intersect(rownames(data), rownames(risk))
  data <- data[sameSample, ]
  risk <- risk[sameSample, ]
  group <- as.vector(risk[,"risk"])
  
  # #PCA分析
  data.class <- rownames(data)
  data.pca <- prcomp(data, scale. = TRUE)
   
  #设置颜色
  color <- ifelse(group == "C1", "#FF9900", "#0066FF")
  
  # Draw PCA plot
  pcaPredict <- predict(data.pca)
  pdf(file = output, width = 7, height = 7)
  par(oma = c(1, 1, 2.5, 1))
  s3d <- scatterplot3d(pcaPredict[, 1:3], pch = 16, color = color, angle = 35)
  legend("top", legend = c("C1", "C2"), pch = 16, inset = -0.2, box.col = "white", xpd = TRUE, horiz = TRUE, col = c("#FF9900", "#0066FF"))
  dev.off()
}

# 所有的PCA分组
#myPCA(input = "normalize.txt", output = "PCA.allGene.pdf")

# Draw PCA plot for stress-related genes
myPCA(input = "内质网应激vennGeneExp.txt", output = "PCA_3D.pdf")

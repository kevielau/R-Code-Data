
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

rm(list=ls()) #一键清除
library(limma)
setwd("E:/BaiduNetdiskDownload/ydh_analysis/AMLreAnalysis/Analysis_Whole/2.dataMerge")        #设置工作目录

#### 1. TCGA 与 GTEx 数据合并 ####

#读取GTEx文件，并整理数据
rt1 <- read.table("symbol_TCGA_normal.txt",sep="\t",header=T,check.names=F)
rt1 <- as.matrix(rt1) #转换成矩阵数据
rownames(rt1) <- rt1[,1] #提取行名
exp1 <- rt1[,2:ncol(rt1)] #提取表达数据
dimnames1 <- list(rownames(exp1),colnames(exp1))
data1 <- matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)
data1 <- avereps(data1) #取均值、去重的操作
data1 <- data1[rowMeans(data1)>0.01,] #将表达均值小于0.01的基因删除

#读取TCGA文件，并整理数据
rt2 <- read.table("symbol_TCGA_tumor.txt",sep="\t",header=T,check.names=F)
rt2 <- as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp2 <- rt2[,2:ncol(rt2)]
dimnames2 <- list(rownames(exp2),colnames(exp2))
data2 <- matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
data2 <- avereps(data2)
data2 <- data2[rowMeans(data2)>0.01,] #将表达均值小于0.01的基因删除

#对基因取交集
sameGene <- intersect( row.names(data1),row.names(data2) ) #取交集
data <- cbind(data1[sameGene,],data2[sameGene,]) #按交集部分合并数据
data[1:4,1:4]

####------------------------------------------------------------------------------------------------------------------
#数据矫正前，需要进行批次效应的验证
group_list=c(rep("control",377),rep("tumor",151))      #设置分组（对照组377例，肿瘤组151例）
table(group_list)

### 1.1 聚类分析
dist_mat <- dist(t(data))
res.hc <- hclust(dist_mat, method = "complete")
#使用dist()函数计算样品之间的距离，并进行聚类分析。

plot(res.hc, labels = colnames(data))
#接着，将聚类结果进行可视化展示。

plot(res.hc, labels = group_list)
#如果显示样品名不够清晰的话，使用分组信息就可以清晰的看到。
#手动保存图片。

### 1.2 PCA分析
library(FactoMineR)
library(factoextra)
data1 <- data #这一步必需，否则后面数据校正会报错！！！
data1 <- PCA(t(data1), graph = FALSE)
fviz_pca_ind(data1,
             geom.ind = "point",
             col.ind = group_list,
             addEllipses = TRUE,
             legend.title = "Groups"
)
ggsave('1.PCA_pre_TCGA.pdf')
#结果显示：两组样品相互混杂，这样会导致差异分析时差异基因数量大大减少。
####------------------------------------------------------------------------------------------------------------------

#数据矫正
outTab <- normalizeBetweenArrays(data)#对两个不同来源的表达数据进行校正和归一化
outTab <- rbind(geneNames=colnames(outTab),outTab)
write.table(outTab,file="1.TCGAmatrix_merge.txt",sep="\t",quote=F,col.names=F)

rm(list=ls()) #一键清除

####=====================================================================================================================
####=====================================================================================================================


#### 2. TARGET 与 GTEx 数据合并 ####

#读取GTEx文件，并整理数据
rt1 <- read.table("symbol_TARGET_normal.txt",sep="\t",header=T,check.names=F)
rt1 <- as.matrix(rt1) #转换成矩阵数据
rownames(rt1) <- rt1[,1] #提取行名
exp1 <- rt1[,2:ncol(rt1)] #提取表达数据
dimnames1 <- list(rownames(exp1),colnames(exp1))
data1 <- matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)
data1 <- avereps(data1) #取均值、去重的操作
data1 <- data1[rowMeans(data1)>0.01,] #将表达均值小于0.01的基因删除

#读取TARGET文件，并整理数据
rt2 <- read.table("symbol_TARGET_tumor.txt",sep="\t",header=T,check.names=F)
rt2 <- as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp2 <- rt2[,2:ncol(rt2)]
dimnames2 <- list(rownames(exp2),colnames(exp2))
data2 <- matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
data2 <- avereps(data2)
data2 <- data2[rowMeans(data2)>0.01,] #将表达均值小于0.01的基因删除

#对基因取交集
sameGene <- intersect( row.names(data1),row.names(data2) ) #取交集
data <- cbind(data1[sameGene,],data2[sameGene,]) #按交集部分合并数据

####------------------------------------------------------------------------------------------------------------------
#数据矫正前，需要进行批次效应的验证
group_list=c(rep("control",755),rep("tumor",358))      #设置分组（对照组377例，肿瘤组151例）
table(group_list)

### 2.1 聚类分析
dist_mat <- dist(t(data))
res.hc <- hclust(dist_mat, method = "complete")
#使用dist()函数计算样品之间的距离，并进行聚类分析。

plot(res.hc, labels = colnames(data))
#接着，将聚类结果进行可视化展示。

plot(res.hc, labels = group_list)
#如果显示样品名不够清晰的话，使用分组信息就可以清晰的看到。
#手动保存图片。

### 2.2 PCA分析
data1 <- data #这一步必需，否则后面数据校正会报错！！！
data1 <- PCA(t(data1), graph = FALSE)
fviz_pca_ind(data1,
             geom.ind = "point",
             col.ind = group_list,
             addEllipses = TRUE,
             legend.title = "Groups"
)
ggsave('1.PCA_pre_TARGET.pdf')
#结果显示：两组样品相互混杂，这样会导致差异分析时差异基因数量大大减少。
####------------------------------------------------------------------------------------------------------------------

#数据矫正
outTab <- normalizeBetweenArrays(data)#对两个不同来源的表达数据进行校正和归一化
outTab <- rbind(geneNames=colnames(outTab),outTab)
write.table(outTab,file="1.TARGETmatrix_merge.txt",sep="\t",quote=F,col.names=F)

rm(list=ls()) #一键清除


#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

rm(list=ls()) #һ�����
library(limma)
setwd("E:/BaiduNetdiskDownload/ydh_analysis/AMLreAnalysis/Analysis_Whole/2.dataMerge")        #���ù���Ŀ¼

#### 1. TCGA �� GTEx ���ݺϲ� ####

#��ȡGTEx�ļ�������������
rt1 <- read.table("symbol_TCGA_normal.txt",sep="\t",header=T,check.names=F)
rt1 <- as.matrix(rt1) #ת���ɾ�������
rownames(rt1) <- rt1[,1] #��ȡ����
exp1 <- rt1[,2:ncol(rt1)] #��ȡ��������
dimnames1 <- list(rownames(exp1),colnames(exp1))
data1 <- matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)
data1 <- avereps(data1) #ȡ��ֵ��ȥ�صĲ���
data1 <- data1[rowMeans(data1)>0.01,] #�������ֵС��0.01�Ļ���ɾ��

#��ȡTCGA�ļ�������������
rt2 <- read.table("symbol_TCGA_tumor.txt",sep="\t",header=T,check.names=F)
rt2 <- as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp2 <- rt2[,2:ncol(rt2)]
dimnames2 <- list(rownames(exp2),colnames(exp2))
data2 <- matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
data2 <- avereps(data2)
data2 <- data2[rowMeans(data2)>0.01,] #�������ֵС��0.01�Ļ���ɾ��

#�Ի���ȡ����
sameGene <- intersect( row.names(data1),row.names(data2) ) #ȡ����
data <- cbind(data1[sameGene,],data2[sameGene,]) #���������ֺϲ�����
data[1:4,1:4]

####------------------------------------------------------------------------------------------------------------------
#���ݽ���ǰ����Ҫ��������ЧӦ����֤
group_list=c(rep("control",377),rep("tumor",151))      #���÷��飨������377����������151����
table(group_list)

### 1.1 �������
dist_mat <- dist(t(data))
res.hc <- hclust(dist_mat, method = "complete")
#ʹ��dist()����������Ʒ֮��ľ��룬�����о��������

plot(res.hc, labels = colnames(data))
#���ţ������������п��ӻ�չʾ��

plot(res.hc, labels = group_list)
#�����ʾ��Ʒ�����������Ļ���ʹ�÷�����Ϣ�Ϳ��������Ŀ�����
#�ֶ�����ͼƬ��

### 1.2 PCA����
library(FactoMineR)
library(factoextra)
data1 <- data #��һ�����裬�����������У���ᱨ��������
data1 <- PCA(t(data1), graph = FALSE)
fviz_pca_ind(data1,
             geom.ind = "point",
             col.ind = group_list,
             addEllipses = TRUE,
             legend.title = "Groups"
)
ggsave('1.PCA_pre_TCGA.pdf')
#�����ʾ��������Ʒ�໥���ӣ������ᵼ�²������ʱ����������������١�
####------------------------------------------------------------------------------------------------------------------

#���ݽ���
outTab <- normalizeBetweenArrays(data)#��������ͬ��Դ�ı������ݽ���У���͹�һ��
outTab <- rbind(geneNames=colnames(outTab),outTab)
write.table(outTab,file="1.TCGAmatrix_merge.txt",sep="\t",quote=F,col.names=F)

rm(list=ls()) #һ�����

####=====================================================================================================================
####=====================================================================================================================


#### 2. TARGET �� GTEx ���ݺϲ� ####

#��ȡGTEx�ļ�������������
rt1 <- read.table("symbol_TARGET_normal.txt",sep="\t",header=T,check.names=F)
rt1 <- as.matrix(rt1) #ת���ɾ�������
rownames(rt1) <- rt1[,1] #��ȡ����
exp1 <- rt1[,2:ncol(rt1)] #��ȡ��������
dimnames1 <- list(rownames(exp1),colnames(exp1))
data1 <- matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)
data1 <- avereps(data1) #ȡ��ֵ��ȥ�صĲ���
data1 <- data1[rowMeans(data1)>0.01,] #�������ֵС��0.01�Ļ���ɾ��

#��ȡTARGET�ļ�������������
rt2 <- read.table("symbol_TARGET_tumor.txt",sep="\t",header=T,check.names=F)
rt2 <- as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp2 <- rt2[,2:ncol(rt2)]
dimnames2 <- list(rownames(exp2),colnames(exp2))
data2 <- matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
data2 <- avereps(data2)
data2 <- data2[rowMeans(data2)>0.01,] #�������ֵС��0.01�Ļ���ɾ��

#�Ի���ȡ����
sameGene <- intersect( row.names(data1),row.names(data2) ) #ȡ����
data <- cbind(data1[sameGene,],data2[sameGene,]) #���������ֺϲ�����

####------------------------------------------------------------------------------------------------------------------
#���ݽ���ǰ����Ҫ��������ЧӦ����֤
group_list=c(rep("control",755),rep("tumor",358))      #���÷��飨������377����������151����
table(group_list)

### 2.1 �������
dist_mat <- dist(t(data))
res.hc <- hclust(dist_mat, method = "complete")
#ʹ��dist()����������Ʒ֮��ľ��룬�����о��������

plot(res.hc, labels = colnames(data))
#���ţ������������п��ӻ�չʾ��

plot(res.hc, labels = group_list)
#�����ʾ��Ʒ�����������Ļ���ʹ�÷�����Ϣ�Ϳ��������Ŀ�����
#�ֶ�����ͼƬ��

### 2.2 PCA����
data1 <- data #��һ�����裬�����������У���ᱨ��������
data1 <- PCA(t(data1), graph = FALSE)
fviz_pca_ind(data1,
             geom.ind = "point",
             col.ind = group_list,
             addEllipses = TRUE,
             legend.title = "Groups"
)
ggsave('1.PCA_pre_TARGET.pdf')
#�����ʾ��������Ʒ�໥���ӣ������ᵼ�²������ʱ����������������١�
####------------------------------------------------------------------------------------------------------------------

#���ݽ���
outTab <- normalizeBetweenArrays(data)#��������ͬ��Դ�ı������ݽ���У���͹�һ��
outTab <- rbind(geneNames=colnames(outTab),outTab)
write.table(outTab,file="1.TARGETmatrix_merge.txt",sep="\t",quote=F,col.names=F)

rm(list=ls()) #һ�����
####=====================================================================================================================####
####=====================================================================================================================####
#### 1. 如何把从 TCGA 下载的 HTseq count 转换成 FPKM ####


#### 1、首先看一下什么是FPKM：
# FPKM（fragments per kilobase million）：fragment per kilobase of transcript per million mapped reads，即每一百万个map上的reads中map到外显子的每一千个碱基上的fragments个数。

# 公式如下：
# FPKM = read counts / (mapped reads(Millions)* exon length) #这里的exon length单位是 kb。

# 这里要先下载的是每个 exon 的长度，为了求长度，我们要先有每个基因组注释文件，可以在这里下载：
# https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files，这是 GDC 官网的网址。进去以后选择这个文件进行下载：
# Annotation Files --- GDC.h38 Flattened GENCODE v22 GFF (used in RNA-Seq alignment and by HTSeq)
#                          genecode.v22.annotation.gtf.gz
#                            md5: 291330bdcff1094bc4d5645de35e0871
#                            file size: 39.0 MB

# 下载后解压备用。


# chooseBioCmirror() #一般选择北京镜像
# chooseCRANmirror() #一般选择兰州镜像
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor") #清华镜像


# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

# BiocManager::install("GenomicFeatures")

# install.packages("tidyverse")



#### 2、计算exon长度，并且保存成dataframe ####

library(GenomicFeatures) #提取exon长度数据
setwd("E:/BaiduNetdiskDownload/R_data_Analysis/2.Papers/Analysis/AML/Analysis_Whole/1.dataPreProcess")
txdb <- makeTxDbFromGFF("gencode.v22.annotation.gtf", format="gtf") #读取上面下载的注释文件

#通过 exonsBy 获取每个基因上的所有外显子的起始位点和终止位点，然后用 reduce 去除掉重叠冗余的部分。
#最后计算长度
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene, function(x){sum(width(reduce(x)))})
#这两个函数太棒了，直接提取注释文件里的基因名（ensemble）和exon长度信息

View(exons_gene_lens) # 得到这个长度后，打开看一下数据集长什么样。
head(exons_gene_lens) # 和 view 函数二选一

#查看有多少基因纳入
length(exons_gene_lens)

#看一下数据类型
class(exons_gene_lens) #查看后发现，是个列表类型的数据：[1] "list"，需要转换成dataframe类型数据

#数据类型转换：list ---> dataframe
exons_gene_lens1 <- as.data.frame(exons_gene_lens)
dim(exons_gene_lens1)#查看转换后的数据维度
str(exons_gene_lens1)#查看数据结构
class(exons_gene_lens1)#查看数据类型

#转换之后，数据列名为基因名，还需要进一步转换：
exons_gene_lens1 <- t(exons_gene_lens1)
dim(exons_gene_lens1)#查看数据维度
head(exons_gene_lens1)#查看数据长什么样，此时行名为基因名。

#删除第一行：
exons_gene_lens2 <- exons_gene_lens1
View(exons_gene_lens2)
exons_gene_lens2 <- as.data.frame(exons_gene_lens2)
rownames(exons_gene_lens2) #查看行名是什么样的
head(exons_gene_lens2)
dim(exons_gene_lens2)


#---------------------------------------------------------------------------------------------------#
#把基因小数点后面的数字去掉：

rownames(exons_gene_lens2) #提取行名
xc <- gsub("\\.(\\.?\\d*)", "", rownames(exons_gene_lens2)) #删除小数点
xc #查看数据
rownames(exons_gene_lens2) = xc #将exon 列表行名等同于 xc，即直接替换exon数据行名后面的小数点。
#---------------------------------------------------------------------------------------------------#


#把 exon 长度表的列名设置为“Length”
colnames(exons_gene_lens2) = "Length"
View(exons_gene_lens2) #查看数据。
write.table(exons_gene_lens2,file = "1.exons_gene_length.txt",sep = "\t",quote = F)
##=======================================================================================================##
#保存文件后，以Excel打开文件，将表头整理一下：第一列加上“id”，“Length”移至第2列，基因名按升序排序后保存。
##=======================================================================================================##

rm(list=ls()) #一键清除

####==========================================================================================================####
####==========================================================================================================####
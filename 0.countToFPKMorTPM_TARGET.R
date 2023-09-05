
setwd("E:/BaiduNetdiskDownload/R_data_Analysis/2.Papers/Analysis/AML/Analysis_Whole/1.dataPreProcess/2.TARGET_cancer")

#### 1.读取表达矩阵 ####
counts <- read.table("mRNAmatrix_TARGET.txt", sep="\t",header=T,check.names=F) #读取表达矩阵
counts[1:4,1:4]
dim(counts)
class(counts) #查看表达矩阵数据类型
rownames(counts) <- counts[,1] #提取行名
exp1 <- counts[,2:ncol(counts)] #提取表达数据

#### 2.读取exon长度数据 ####
exon_gene_length <- read.table("1.exons_gene_length.txt",sep = "\t",header = T,check.names = F) #读取exon长度文件
exon_gene_length[1:4,1:2] #exon长度表只有两列。
dim(exon_gene_length)
class(exon_gene_length) #查看exon长度表的数据类型
rownames(exon_gene_length) <- exon_gene_length[,1] #提取行名
exp2 <- exon_gene_length[,2:ncol(exon_gene_length)] #提取表达数据

# 去除exon文件的基因名中的小数点，方便后面与基因表达矩阵合并
xc <- gsub("\\.(\\.?\\d*)", "", rownames(exon_gene_length)) #删除小数点,记住该函数！！！！！！！！！！！
xc #查看数据
rownames(exon_gene_length) = xc #将exon 列表行名等同于 xc，即直接替换exon数据行名后面的小数点。

#------------------------------------------------------------------------------------#
#可以看出表达矩阵和exon长度数据的基因名的区别：前者ensemble名没有小数点，后者有。
#所以为了合并，exon 的基因名也需要处理一下：去除小数点。
#class()函数先检查两种数据是否都是data.frame，如果不是，则需做如下转换：
               #counts <- as.data.frame(counts)
               #exon_gene_length <- as.data.frame(exon_gene_length)
#------------------------------------------------------------------------------------#


#### 3.对基因取交集 ####
sameGene <- intersect( row.names(counts),row.names(exon_gene_length) ) #取交集
data <- cbind(counts[sameGene,],exon_gene_length[sameGene,]) #按交集部分合并数据
dim(data)
class(data)

#### 4.数据保存 ####
outTab <- rbind(geneNames=colnames(data),data)
write.table(outTab,file="2.count_with_length.txt",sep="\t",quote=F,col.names=F)
#写入文件后，用Excel打开文件，删除多余的列。

##=========================================================================================================================##
##=========================================================================================================================##
##=========================================================================================================================##

#### 5.计算 fpkm 值：

# 先把 length 除以 1000，就是上面公式里说的单位要 kb
count_with_length <- read.table("2.count_with_length.txt",sep="\t",header=T,check.names=F,row.names = 1) 
#row.names = 1定义行名之后，第1-358列就是样本名，第359列是Length名。
#如果没定义行名，则2-359列为样本名
dim(count_with_length)
class(count_with_length) #data.frame类型数据
kb <- count_with_length$Length/1000
kb #查看数据

# 把合并矩阵里的前 151 列的数值都除以 kb：这个根据自己的数据来，本示例数据样本数为151，即1:151。
countdata <- count_with_length[,1:(ncol(count_with_length)-1)] #倒数第一列是exon的长度数据，故用总列数-1
rpk <- countdata/kb
rpk #查看数据

#计算 FPKM 值
fpkmMatrix <- t(t(rpk)/colSums(countdata)*10^6)
head(fpkmMatrix) #查看数据
fpkmMatrix[1:4,1:4]

#保存 FPKM 矩阵
write.table(fpkmMatrix, file="3.fpkmMatrix.txt",sep="\t",quote=F)
#保存后需要修改表头。

####===========================================================================================================####
####===========================================================================================================####


## 二、接下来我们要进行FPKM转为TPM的转换

# 1. 其实原理性的东西大家大可不必去深究，反正就记住TPM的数据要比FPKM更准，能使用TPM还是要尽量使用TPM，转换公式如下：

fpkmToTpm <- function(fpkm)
{ 
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


#数据类型的固定转换公式如下，按需索取：
#参考网站：https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/

#countToTpm <- function(counts, effLen)                 ##count值转换为TPM值
#{
# rate <- log(counts) - log(effLen)
#denom <- log(sum(exp(rate)))
#exp(rate - denom + log(1e6))
#}

#countToFpkm <- function(counts, effLen)                ##count值转换为FPKM值
#{
# N <- sum(counts)
#exp( log(counts) + log(1e9) - log(effLen) - log(N) )
#}

#fpkmToTpm <- function(fpkm)                            ##FPKM值转换为TPM值
#{
# exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
#}

#countToEffCounts <- function(counts, len, effLen)      ##count值转换为effect count值
#{
# counts * (len / effLen)
#}

# 2. 接下来我们使用apply函数对FPKM数据的按照列进行转换，一定要注意这里是按照列进行转换！！！

tpmMatrix <- apply(fpkmMatrix, 2, fpkmToTpm) ##这里是按列名，即样本名，来进行转换的。
#其他数据类型转换，也用apply函数进行。

write.table(tpmMatrix, file="4.tpmMatrix.txt",sep="\t",quote=F)

rm(list=ls()) #一键清除

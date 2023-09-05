
#### 1. TCGA 数据 ####

#### 1.1 PCA 分析 ####
rm(list = ls())
options(stringsAsFactors = F)

setwd("E:/BaiduNetdiskDownload/ydh_analysis/AMLreAnalysis/Analysis_Whole/2.dataMerge")
rt=read.table("1.TCGAmatrix_merge.txt",sep="\t",header=T,check.names=F,row.names = 1)  #读取合并后的表达数据
#该数据已经过去重处理，所以不存在重复基因名，可直接定义行名（row.names = 1）

group_list=c(rep("control",377),rep("tumor",151))      #设置分组（对照组377例，肿瘤组151例）

table(group_list) #统计肿瘤和对照数量
dat=rt
# 每次都要检测数据
dat[1:4,1:4]
## 下面是画PCA的必须操作，需要看说明书。
dat=t(dat)                  #矩阵转置
dat=as.data.frame(dat)      #转换成数据框
dat=cbind(dat,group_list)   #合并分组
library("FactoMineR")
library("factoextra") 
# The variable group_list (index = ncol(dat)) is removed
# before PCA analysis
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE) #PCA分析
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list, # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('2.PCA_post_TCGA.pdf')

#### 1.2 heatmap绘制 ####

dat=log(rt+1)   #######################数据的 log2 转换。
# 每次都要检测数据
dat[1:4,1:4]
cg=names(tail(sort(apply(dat,1,sd)),1000))# 选择SD值最大的前1000个基因进行热图绘制 
mat=dat[cg,]
library(pheatmap)
pheatmap(mat,show_colnames =F,show_rownames = F)
n=t(scale(t(mat)))
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(g=group_list)
rownames(ac)=colnames(n)
#可以看到，对照与肿瘤之间的表达量可以很好地区分开。
pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac,filename = '2.heatmap_post_top1000SD_TCGA.pdf')


####-------------------------------------------------------------------------------------------------------

#### 2.TARGET数据 ####
rm(list = ls())
options(stringsAsFactors = F)

rt=read.table("1.TARGETmatrix_merge.txt",sep="\t",header=T,check.names=F,row.names = 1)  #读取合并后的表达数据
#该数据已经过去重处理，所以不存在重复基因名，可直接定义行名（row.names = 1）

group_list=c(rep("control",755),rep("tumor",358))      #设置分组（对照组377例，肿瘤组151例）

table(group_list) #统计肿瘤和对照数量
dat=rt #建立一个相同的新数据集，方便后续继续分析。
# 每次都要检测数据
dat[1:4,1:4]
## 下面是画PCA的必须操作，需要看说明书。
dat=t(dat)                  #矩阵转置
dat=as.data.frame(dat)      #转换成数据框
dat=cbind(dat,group_list)   #合并分组
library("FactoMineR")
library("factoextra") 
# The variable group_list (index = ncol(dat)) is removed
# before PCA analysis
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE) #PCA分析
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list, # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('2.PCA_post_TARGET.pdf')

#### 1.2 heatmap绘制 ####

dat=log(rt+1)   #######################数据的 log2 转换。
# 每次都要检测数据
dat[1:4,1:4]
cg=names(tail(sort(apply(dat,1,sd)),1000))# 选择SD值最大的前1000个基因进行热图绘制 
mat=dat[cg,]
library(pheatmap)
pheatmap(mat,show_colnames =F,show_rownames = F)
n=t(scale(t(mat)))
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(g=group_list)
rownames(ac)=colnames(n)
#可以看到，对照与肿瘤之间的表达量可以很好地区分开。
pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac,filename = '2.heatmap_post_top1000SD_TARGET.pdf')

rm(list = ls())

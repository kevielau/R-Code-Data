
#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE", version = "3.8")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("clusterProfiler", version = "3.8")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("enrichplot", version = "3.8")


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

setwd("E:/BaiduNetdiskDownload/ydh_analysis/AML/Analysis_Whole/5.TCGA_TME/19.KEGG")
rt=read.table("id_TCGA.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,] #过滤NA值
gene=rt$entrezID

#kegg富集分析
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG_TCGA.txt",sep="\t",quote=F,row.names = F)

#柱状图
tiff(file="barplot_TCGA.tiff",width = 20,height = 12,units ="cm",compression="lzw",bg="white",res=600)
barplot(kk, drop = TRUE, showCategory = 30)
dev.off()

#气泡图
tiff(file="dotplot_TCGA.tiff",width = 20,height = 12,units ="cm",compression="lzw",bg="white",res=600)
dotplot(kk, showCategory = 30)
dev.off()

rm(list = ls()) #一键清除


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

setwd("E:/BaiduNetdiskDownload/ydh_analysis/AML/Analysis_Whole/5.TCGA_TME/18.GO")
rt=read.table("id_TCGA.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]
gene=rt$entrezID

#GO富集分析
kk <- enrichGO(gene = gene, #按照ensembl名进行富集分析
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all", #是对bp/cc/mf全部做注释
               readable =T)
write.table(kk,file="GO_TCGA.txt",sep="\t",quote=F,row.names = F)

#柱状图
tiff(file="barplot_TCGA.tiff",width = 26,height = 20,units ="cm",compression="lzw",bg="white",res=600)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#气泡图
tiff(file="dotplot_TCGA.tiff",width = 26,height = 20,units ="cm",compression="lzw",bg="white",res=600)
dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

rm(list = ls()) #一键清除

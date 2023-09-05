
#install.packages("pheatmap")

setwd("E:/BaiduNetdiskDownload/TARGET_database/TARGET/09.pheatmap")      #设置工作目录
rt=read.table("diffGeneExp.txt",sep="\t",header=T,row.names=1,check.names=F) #用的是差异表达基因的表达数据
rt=log2(rt+0.001) #log2转置

library(pheatmap)
Type=c(rep("Female",172),rep("Male",186))    #修改女性和男性样品数目
names(Type)=colnames(rt)
Type=as.data.frame(Type)

tiff(file="heatmap_ydh.tiff",
       width = 70,            #图片的宽度
       height =40,            #图片的高度
       units ="cm",
       compression="lzw",
       bg="white",
       res=600) #分辨率
pheatmap(rt, 
         annotation=Type, #按分组注释
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F, #不聚类，这样分组就不会被打乱顺序
#         scale="row",
         fontsize_row=7,
         fontsize_col=5)
dev.off()

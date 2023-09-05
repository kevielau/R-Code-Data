
#install.packages("pheatmap")

setwd("E:/BaiduNetdiskDownload/TARGET_database/TARGET/17.pheatmap")
rt=read.table("risk.txt",sep="\t",header=T,row.names=1,check.names=F)
rt=rt[order(rt$riskScore),]#按风险值排序
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1) #转置

rt1=log2(rt1+1)
library(pheatmap)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt) #注释高、低风险

tiff(file="heatmap_ydh.tiff",width = 50,height = 20,units ="cm",compression="lzw",bg="white",res=300)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE, #不聚类
         fontsize_row=11,
         fontsize_col=6,
         color = colorRampPalette(c("green", "black", "red"))(50) )
dev.off()

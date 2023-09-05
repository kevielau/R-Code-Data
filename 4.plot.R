
#install.packages("pheatmap")

setwd("E:/BaiduNetdiskDownload/TARGET_database/TARGET/09.pheatmap")      #���ù���Ŀ¼
rt=read.table("diffGeneExp.txt",sep="\t",header=T,row.names=1,check.names=F) #�õ��ǲ���������ı�������
rt=log2(rt+0.001) #log2ת��

library(pheatmap)
Type=c(rep("Female",172),rep("Male",186))    #�޸�Ů�Ժ�������Ʒ��Ŀ
names(Type)=colnames(rt)
Type=as.data.frame(Type)

tiff(file="heatmap_ydh.tiff",
       width = 70,            #ͼƬ�Ŀ���
       height =40,            #ͼƬ�ĸ߶�
       units ="cm",
       compression="lzw",
       bg="white",
       res=600) #�ֱ���
pheatmap(rt, 
         annotation=Type, #������ע��
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F, #�����࣬��������Ͳ��ᱻ����˳��
#         scale="row",
         fontsize_row=7,
         fontsize_col=5)
dev.off()
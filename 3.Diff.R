
options(stringsAsFactors = F)

library("limma")
setwd("E:/BaiduNetdiskDownload/ydh_analysis/AMLreAnalysis/Analysis_Whole/2_1limmaDiff")       #设置工作目录

##========================================================================================================================#

#### 1.TCGA 差异分析 ####

inputFile="1.TCGAmatrix_merge.txt"                                        #输入文件
fdrFilter=0.01                                                    #fdr临界值
logFCfilter=2                                                     #logFC临界值
conNum=377                                                        #normal组样品数目
treatNum=151                                                      #tumor组样品数目

#读取输入文件
outTab=data.frame() #创建数据框
grade=c(rep(1,conNum),rep(2,treatNum))
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data) #上述代码是去重操作
data=data[rowMeans(data)>0,]

#差异分析
for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  #差异分析的检验方法
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr") #计算FDR值
outTab=cbind(outTab,fdr=fdr)

#输出所有基因的差异情况
write.table(outTab,file="all_TCGA.xls",sep="\t",row.names=F,quote=F)

#输出差异表格
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff_TCGA.xls",sep="\t",row.names=F,quote=F)

#绘制热图需要的文件
heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="diffExp_TCGA.txt",sep="\t",col.names=F,quote=F)

#绘制火山图
pdf(file="vol_TCGA.pdf",height=8,width=10)
xMax=max(abs(as.numeric(as.vector(outTab$logFC))))
yMax=max(-log10(outTab$fdr))+1
plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="red",cex=1)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="green",cex=1)
abline(v=0,lty=2,lwd=3)
dev.off()

#绘制差异基因热图
library(pheatmap)
hmExp=data[as.vector(outDiff[,1]),]
hmExp=log2(hmExp+0.01)
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap_TCGA.pdf",height=8,width=10)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(50),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = T,
         fontsize = 12,
         fontsize_row=3,
         fontsize_col=3)
dev.off()

rm(list = ls()) #一键清除


#===================================================================================================================================#
#===================================================================================================================================#

#### 2.TARGET 差异分析 ####

inputFile="1.TARGETmatrix_merge.txt"                                      #输入文件
fdrFilter=0.01                                                    #fdr临界值
logFCfilter=2                                                     #logFC临界值
conNum=755                                                        #normal组样品数目
treatNum=358                                                      #tumor组样品数目

#读取输入文件
outTab=data.frame() #创建数据框
grade=c(rep(1,conNum),rep(2,treatNum))
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data) #上述代码是去重操作
data=data[rowMeans(data)>0,]

#差异分析
for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  #差异分析的检验方法
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr") #计算FDR值
outTab=cbind(outTab,fdr=fdr)

#输出所有基因的差异情况
write.table(outTab,file="all_TARGET.xls",sep="\t",row.names=F,quote=F)

#输出差异表格
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff_TARGET.xls",sep="\t",row.names=F,quote=F)

#绘制热图需要的文件
heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="diffExp_TARGET.txt",sep="\t",col.names=F,quote=F)

#绘制火山图
pdf(file="vol_TARGET.pdf",height=8,width=10)
xMax=max(abs(as.numeric(as.vector(outTab$logFC))))
yMax=max(-log10(outTab$fdr))+1
plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="red",cex=1)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="green",cex=1)
abline(v=0,lty=2,lwd=3)
dev.off()

#绘制差异基因热图
library(pheatmap)
hmExp=data[as.vector(outDiff[,1]),]
hmExp=log2(hmExp+0.01)
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap_TARGET.pdf",height=8,width=10)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(50),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = T,
         fontsize = 12,
         fontsize_row=3,
         fontsize_col=3)
dev.off()

rm(list = ls()) #一键清除


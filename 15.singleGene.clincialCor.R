#### 1. TCGA ���� ####

#### 1.1 VDR ���� ####
inputFile="singleGeneClinical_VDR_TCGA.txt"                                     #�����ļ�
setwd("E:/BaiduNetdiskDownload/ydh_analysis/AML/Analysis_Whole/6_singleGene_revised/14.clinicalCor")         #�޸Ĺ���Ŀ¼
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
clinical="stage"                                                       #�����ٴ�����
gene="VDR"                                                            #�����������

xlabel=vector()
tab1=table(rt[,clinical])
labelNum=length(tab1)
dotCol=c(2,3)  #�������ӵ���ɫ
if(labelNum==3){
	dotCol=c(2,3,4)
}
if(labelNum==4){
	dotCol=c(2,3,4,5)
}
if(labelNum>4){
	dotCol=rainbow(labelNum)
}
for(i in 1:labelNum){
  xlabel=c(xlabel,names(tab1[i]) ) #����stage������
}

outTab=data.frame()

rt1=rbind(expression=rt[,gene],clinical=rt[,clinical])
rt1=as.matrix(t(rt1))
if(labelNum==2){
    wilcoxTest<-wilcox.test(expression ~ clinical, data=rt1) #�����Ƚ���wilcox����
}else{
    wilcoxTest<-kruskal.test(expression ~ clinical, data = rt1)}#�����Ƚ���kruskal����
pValue=wilcoxTest$p.value
outTab=rbind(outTab,cbind(gene=gene,pVal=pValue))
pval=0
if(pValue<0.001){
    pval=signif(pValue,4)
    pval=format(pval, scientific = TRUE)
}else{
    pval=round(pValue,3)
}
  
b = boxplot(expression ~ clinical, data = rt1,outline = FALSE, plot=F) 
yMin=min(b$stats)
yMax = max(b$stats/5+b$stats)
ySeg = max(b$stats/10+b$stats)
ySeg2 = max(b$stats/12+b$stats)
n = ncol(b$stats)

outPdf=paste0(gene,".",clinical,".pdf")
pdf(file=outPdf,
    width=9,
    height=6,)
par(mar = c(4,7,3,3))
boxplot(expression ~ clinical, data = rt1,names=xlabel,
    ylab = paste0(gene," expression"),col=dotCol,
    cex.main=1.6, cex.lab=1.4, cex.axis=1.3,ylim=c(yMin,yMax),outline = FALSE)
segments(1,ySeg, n,ySeg);
segments(1,ySeg, 1,ySeg2)
segments(n,ySeg, n,ySeg2)
text((1+n)/2,ySeg,labels=paste0("p=",pval),cex=1.5,pos=3)
dev.off()

rm(list = ls())


#### 1.2 S100A6 ���� ####
inputFile="singleGeneClinical_S100A6_TCGA.txt"                                     #�����ļ�
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
clinical="stage"                                                       #�����ٴ�����
gene="S100A6"                                                            #�����������

xlabel=vector()
tab1=table(rt[,clinical])
labelNum=length(tab1)
dotCol=c(2,3)  #�������ӵ���ɫ
if(labelNum==3){
  dotCol=c(2,3,4)
}
if(labelNum==4){
  dotCol=c(2,3,4,5)
}
if(labelNum>4){
  dotCol=rainbow(labelNum)
}
for(i in 1:labelNum){
  xlabel=c(xlabel,names(tab1[i]) ) #����stage������
}

outTab=data.frame()

rt1=rbind(expression=rt[,gene],clinical=rt[,clinical])
rt1=as.matrix(t(rt1))
if(labelNum==2){
  wilcoxTest<-wilcox.test(expression ~ clinical, data=rt1) #�����Ƚ���wilcox����
}else{
  wilcoxTest<-kruskal.test(expression ~ clinical, data = rt1)}#�����Ƚ���kruskal����
pValue=wilcoxTest$p.value
outTab=rbind(outTab,cbind(gene=gene,pVal=pValue))
pval=0
if(pValue<0.001){
  pval=signif(pValue,4)
  pval=format(pval, scientific = TRUE)
}else{
  pval=round(pValue,3)
}

b = boxplot(expression ~ clinical, data = rt1,outline = FALSE, plot=F) 
yMin=min(b$stats)
yMax = max(b$stats/5+b$stats)
ySeg = max(b$stats/10+b$stats)
ySeg2 = max(b$stats/12+b$stats)
n = ncol(b$stats)

outPdf=paste0(gene,".",clinical,".pdf")
pdf(file=outPdf,
    width=9,
    height=6,)
par(mar = c(4,7,3,3))
boxplot(expression ~ clinical, data = rt1,names=xlabel,
        ylab = paste0(gene," expression"),col=dotCol,
        cex.main=1.6, cex.lab=1.4, cex.axis=1.3,ylim=c(yMin,yMax),outline = FALSE)
segments(1,ySeg, n,ySeg);
segments(1,ySeg, 1,ySeg2)
segments(n,ySeg, n,ySeg2)
text((1+n)/2,ySeg,labels=paste0("p=",pval),cex=1.5,pos=3)
dev.off()

rm(list = ls())


####--------------------------------------------------------------------------------------

#### 2. TARGET ���� ####

#### 2.1 VDR ���� ####
inputFile="singleGeneClinical_VDR_TARGET.txt"                                     #�����ļ�
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
clinical="stage"                                                       #�����ٴ�����
gene="VDR"                                                            #�����������

xlabel=vector()
tab1=table(rt[,clinical])
labelNum=length(tab1)
dotCol=c(2,3)  #�������ӵ���ɫ
if(labelNum==3){
  dotCol=c(2,3,4)
}
if(labelNum==4){
  dotCol=c(2,3,4,5)
}
if(labelNum>4){
  dotCol=rainbow(labelNum)
}
for(i in 1:labelNum){
  xlabel=c(xlabel,names(tab1[i]) ) #����stage������
}

outTab=data.frame()

rt1=rbind(expression=rt[,gene],clinical=rt[,clinical])
rt1=as.matrix(t(rt1))
if(labelNum==2){
  wilcoxTest<-wilcox.test(expression ~ clinical, data=rt1) #�����Ƚ���wilcox����
}else{
  wilcoxTest<-kruskal.test(expression ~ clinical, data = rt1)}#�����Ƚ���kruskal����
pValue=wilcoxTest$p.value
outTab=rbind(outTab,cbind(gene=gene,pVal=pValue))
pval=0
if(pValue<0.001){
  pval=signif(pValue,4)
  pval=format(pval, scientific = TRUE)
}else{
  pval=round(pValue,3)
}

b = boxplot(expression ~ clinical, data = rt1,outline = FALSE, plot=F) 
yMin=min(b$stats)
yMax = max(b$stats/5+b$stats)
ySeg = max(b$stats/10+b$stats)
ySeg2 = max(b$stats/12+b$stats)
n = ncol(b$stats)

outPdf=paste0(gene,".",clinical,".pdf")
pdf(file=outPdf,
    width=9,
    height=6,)
par(mar = c(4,7,3,3))
boxplot(expression ~ clinical, data = rt1,names=xlabel,
        ylab = paste0(gene," expression"),col=dotCol,
        cex.main=1.6, cex.lab=1.4, cex.axis=1.3,ylim=c(yMin,yMax),outline = FALSE)
segments(1,ySeg, n,ySeg);
segments(1,ySeg, 1,ySeg2)
segments(n,ySeg, n,ySeg2)
text((1+n)/2,ySeg,labels=paste0("p=",pval),cex=1.5,pos=3)
dev.off()

rm(list = ls())


#### 2.2 S100A6 ���� ####
inputFile="singleGeneClinical_S100A6_TARGET.txt"                                     #�����ļ�
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
clinical="stage"                                                       #�����ٴ�����
gene="S100A6"                                                            #�����������

xlabel=vector()
tab1=table(rt[,clinical])
labelNum=length(tab1)
dotCol=c(2,3)  #�������ӵ���ɫ
if(labelNum==3){
  dotCol=c(2,3,4)
}
if(labelNum==4){
  dotCol=c(2,3,4,5)
}
if(labelNum>4){
  dotCol=rainbow(labelNum)
}
for(i in 1:labelNum){
  xlabel=c(xlabel,names(tab1[i]) ) #����stage������
}

outTab=data.frame()

rt1=rbind(expression=rt[,gene],clinical=rt[,clinical])
rt1=as.matrix(t(rt1))
if(labelNum==2){
  wilcoxTest<-wilcox.test(expression ~ clinical, data=rt1) #�����Ƚ���wilcox����
}else{
  wilcoxTest<-kruskal.test(expression ~ clinical, data = rt1)}#�����Ƚ���kruskal����
pValue=wilcoxTest$p.value
outTab=rbind(outTab,cbind(gene=gene,pVal=pValue))
pval=0
if(pValue<0.001){
  pval=signif(pValue,4)
  pval=format(pval, scientific = TRUE)
}else{
  pval=round(pValue,3)
}

b = boxplot(expression ~ clinical, data = rt1,outline = FALSE, plot=F) 
yMin=min(b$stats)
yMax = max(b$stats/5+b$stats)
ySeg = max(b$stats/10+b$stats)
ySeg2 = max(b$stats/12+b$stats)
n = ncol(b$stats)

outPdf=paste0(gene,".",clinical,".pdf")
pdf(file=outPdf,
    width=9,
    height=6,)
par(mar = c(4,7,3,3))
boxplot(expression ~ clinical, data = rt1,names=xlabel,
        ylab = paste0(gene," expression"),col=dotCol,
        cex.main=1.6, cex.lab=1.4, cex.axis=1.3,ylim=c(yMin,yMax),outline = FALSE)
segments(1,ySeg, n,ySeg);
segments(1,ySeg, 1,ySeg2)
segments(n,ySeg, n,ySeg2)
text((1+n)/2,ySeg,labels=paste0("p=",pval),cex=1.5,pos=3)
dev.off()

rm(list = ls())

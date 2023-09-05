
#install.packages("glmnet")
#install.packages("survival")


library("glmnet")
library("survival")

setwd("E:/BaiduNetdiskDownload/ydh_analysis/AMLreAnalysis/Analysis_Whole/4.immuneGene/11.lasso_multicox/06.lasso")                #设置工作目录

#### 1. TARGET 数据 ####

rt=read.table("uniSigExp_TARGET.txt",header=T,sep="\t",row.names=1,check.names=F)     #读取文件
rt$futime[rt$futime<=0]=1 #如果随访时间有0的，lasso分析会报错。所以这里把小于等于0的值设为1.

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

fit <- glmnet(x, y, family = "cox", maxit = 1000)
#lasso分析。maxit = 1000随机模拟1000次。
pdf("lambda_TARGET.pdf")
plot(fit, xvar = "lambda", label = TRUE) #"lambda"是惩罚系数，越大，标准越严格。
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("cvfit_TARGET.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")#标出线条
dev.off()

#提取误差最小、系数不为0的基因提取出来
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lassoSigExp_TARGET.txt",sep="\t",row.names=F,quote=F)#输出文件

rm(list = ls())#一键清除

####================================================================================================

#### 2. TCGA 数据 ####

rt=read.table("uniSigExp_TCGA.txt",header=T,sep="\t",row.names=1,check.names=F)     #读取文件
rt$futime[rt$futime<=0]=1 #如果随访时间有0的，lasso分析会报错。所以这里把小于等于0的值设为1.

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

fit <- glmnet(x, y, family = "cox", maxit = 1000)
#lasso分析。maxit = 1000随机模拟1000次。
pdf("lambda_TCGA.pdf")
plot(fit, xvar = "lambda", label = TRUE) #"lambda"是惩罚系数，越大，标准越严格。
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("cvfit_TCGA.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")#标出线条
dev.off()

#提取误差最小、系数不为0的基因提取出来
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lassoSigExp_TCGA.txt",sep="\t",row.names=F,quote=F)#输出文件

rm(list = ls())#一键清除


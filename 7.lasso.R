
#install.packages("glmnet")
#install.packages("survival")


library("glmnet")
library("survival")

setwd("E:/BaiduNetdiskDownload/ydh_analysis/AMLreAnalysis/Analysis_Whole/4.immuneGene/11.lasso_multicox/06.lasso")                #���ù���Ŀ¼

#### 1. TARGET ���� ####

rt=read.table("uniSigExp_TARGET.txt",header=T,sep="\t",row.names=1,check.names=F)     #��ȡ�ļ�
rt$futime[rt$futime<=0]=1 #������ʱ����0�ģ�lasso�����ᱨ�������������С�ڵ���0��ֵ��Ϊ1.

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

fit <- glmnet(x, y, family = "cox", maxit = 1000)
#lasso������maxit = 1000���ģ��1000�Ρ�
pdf("lambda_TARGET.pdf")
plot(fit, xvar = "lambda", label = TRUE) #"lambda"�ǳͷ�ϵ����Խ�󣬱�׼Խ�ϸ�
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("cvfit_TARGET.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")#�������
dev.off()

#��ȡ�����С��ϵ����Ϊ0�Ļ�����ȡ����
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lassoSigExp_TARGET.txt",sep="\t",row.names=F,quote=F)#����ļ�

rm(list = ls())#һ�����

####================================================================================================

#### 2. TCGA ���� ####

rt=read.table("uniSigExp_TCGA.txt",header=T,sep="\t",row.names=1,check.names=F)     #��ȡ�ļ�
rt$futime[rt$futime<=0]=1 #������ʱ����0�ģ�lasso�����ᱨ�������������С�ڵ���0��ֵ��Ϊ1.

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

fit <- glmnet(x, y, family = "cox", maxit = 1000)
#lasso������maxit = 1000���ģ��1000�Ρ�
pdf("lambda_TCGA.pdf")
plot(fit, xvar = "lambda", label = TRUE) #"lambda"�ǳͷ�ϵ����Խ�󣬱�׼Խ�ϸ�
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("cvfit_TCGA.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")#�������
dev.off()

#��ȡ�����С��ϵ����Ϊ0�Ļ�����ȡ����
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lassoSigExp_TCGA.txt",sep="\t",row.names=F,quote=F)#����ļ�

rm(list = ls())#һ�����

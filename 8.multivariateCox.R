
#install.packages("survival")

library(survival)
setwd("E:/BaiduNetdiskDownload/TARGET_database/TARGET/14.multiCox")
rt=read.table("multiInput.txt",header=T,sep="\t",check.names=F,row.names=1)
rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1) #第3列开始为基因名，对基因表达量进行log2处理
rt[,"futime"]=rt[,"futime"]/365
cox <- coxph(Surv(futime, fustat) ~ ., data = rt)
cox=step(cox,direction = "both")
riskScore=predict(cox,type="risk",newdata=rt)
summary=summary(cox)
coxGene=rownames(summary$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
    file="risk_ydh.txt",
    sep="\t",
    quote=F,
    row.names=F)
write.table(cbind(id=coxGene,summary$coefficients),
    file="coxResult_ydh.xls",
    sep="\t",
    quote=F,
    row.names=F)
cox

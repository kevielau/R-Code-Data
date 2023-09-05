
setwd("E:/BaiduNetdiskDownload/TARGET_database/TARGET/12.univariateCox")
outTab=data.frame()

library(survival)
rt=read.table("clinicalExp.txt",header=T,sep="\t",row.names=1,check.names=F)
rt1=log2(rt[,3:ncol(rt)]+1)
rt=cbind(rt[,1:2],rt1)

for(i in colnames(rt[,3:ncol(rt)])){  #第3列开始，是基因名。是对基因表达进行cox分析。
 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)  
 #生存时间和状态为y值，基因表达（rt[,i]）为x值。
 #每次拿出一个基因进行分析，一直循环到最后一个基因。得到的是单因素cox。
 coxSummary = summary(cox) #保留cox分析数据
 outTab=rbind(outTab,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"], #没分析一个基因，就保存到“outTab”中。
 z=coxSummary$coefficients[,"z"],
 pvalue=coxSummary$coefficients[,"Pr(>|z|)"])) #p值是通过z值计算出来的。
}

write.table(outTab,file="univariateCox_ydh.xls",sep="\t",row.names=F,quote=F)

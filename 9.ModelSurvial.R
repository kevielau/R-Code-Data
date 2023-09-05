
#install.packages("survival")

setwd("E:/BaiduNetdiskDownload/TARGET_database/TARGET/15.survival")   
library(survival)
rt=read.table("risk.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt) #高低风险两组比较。
pValue=1-pchisq(diff$chisq,df=1)
#pValue=round(pValue,3)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE) #科学计数法

fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)    #查看5年生存率
tiff(file="survival_ydh.tiff",
       width = 14,            
       height =14,            
       units ="cm",
       compression="lzw",
       bg="white",
       res=600)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""),
     mark.time=T)
legend("topright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()

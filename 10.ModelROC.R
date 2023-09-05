
#install.packages("survivalROC")

library(survivalROC)
setwd("E:/BaiduNetdiskDownload/TARGET_database/TARGET/16.ROC")
rt=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1) #风险值文件
#check.names=F，不能检查表头，因为样本名含有“-”符号，检查的话会被软件自动改为“.”
tiff(file="ROC_ydh.tiff",
       width = 14,            
       height =14,            
       units ="cm",
       compression="lzw",
       bg="white",
       res=600)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                #根据风险评分来标记。
      predict.time =5, method="KM") #5年预测时间（根据需要，可改为1或3年）。
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
  xlab="False positive rate", ylab="True positive rate",
  main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()



setwd("E:/BaiduNetdiskDownload/TARGET_database/TARGET/18.riskScore")

rt=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"] #按 "risk" 列分类
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10 #风险值大于10时，定义为10，以防整体趋势偏差太大。
tiff(file="riskScore_ydh.tiff",width = 20, height = 12,units ="cm",compression="lzw",bg="white",res=600)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
     rep("red",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2) #在高低风险临界值处画横线
dev.off()

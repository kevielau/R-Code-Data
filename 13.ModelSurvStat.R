
setwd("E:/BaiduNetdiskDownload/TARGET_database/TARGET/19.survStat")

rt=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
color=as.vector(rt$fustat)
color[color==1]="red" #死亡病例用红色表示
color[color==0]="green" #活着的病人绿色表示
tiff(file="survStat_ydh.tiff",width = 20, height = 12,units ="cm",compression="lzw",bg="white",res=600)
plot(rt$futime,
     pch=19, #画点
     xlab="Patients (increasing risk socre)", #按打分来排序。
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
dev.off()

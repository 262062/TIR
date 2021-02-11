
library(survival)

rt=read.table("Time.txt",header=T,sep="\t",check.names=F)   
rt$futime=rt$futime/365                                           
rt$cluster=factor(rt$cluster, levels=c("L","M","H"))

diff=survdiff(Surv(futime, fustat) ~cluster,data = rt)
pValue=1-pchisq(diff$chisq,df=2)
if(pValue<0.001){
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
}else{
  pValue=round(pValue,3)
}

fit <- survfit(Surv(futime, fustat) ~ cluster, data = rt)

pdf(file="survival.pdf",width = 5.5,height =5)
plot(fit, 
     lwd=2,
     col=c("green","blue","red"),
     xlab="Time (year)",
     mark.time=T,
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""))
legend("topright", 
     c("L","M","H"), 
     lwd=2, 
     col=c("green","blue","red"))

summary(fit)


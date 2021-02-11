
library("glmnet")
library("survival")
rt=read.table("tcgaSig.txt",header=T,sep="\t",row.names=1,check.names=F)   
rt$futime[rt$futime<=0]=1
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("lam.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("fit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lasso=row.names(coef)[index]
lasso=c("futime","fustat",lassoGene)
lassoSig=rt[,lassoGene]
lassoSig=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSig,file="lassoSig.txt",sep="\t",row.names=F,quote=F)

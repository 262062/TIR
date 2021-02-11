
library(ggpubr)
rt=read.table("PDL1CTLA4exp.txt",sep="\t",header=T,row.names=1,check.names=F)   
Type=read.table("cluster.Immunity.txt",sep="\t",check.names=F,row.names=1,header=F)
Type=Type[order(Type[,2]),]
rt=t(rt[,row.names(Type)])
data=data.frame()
for(i in colnames(rt)){
  data=rbind(data,cbind(expression=log2(rt[,i]+1),gene=i,Subtype=as.vector(Type[,2])))
}
write.table(data,file="data.txt",sep="\t",row.names=F,quote=F)
data=read.table("MHCIdata.txt",sep="\t",header=T,check.names=F)     
data$Subtype=factor(data$Subtype, levels=c("Immunity_L","Immunity_M","Immunity_H"))
p=ggboxplot(data, x="gene", y="expression", color = "Subtype", 
     ylab="Gene expression (log2(FPKM+1))",
     xlab="",
     palette = c("chartreuse4","blue","red") )
}
#p=p+rotate_x_text(100)
axis.title.x=element_text(size=80)
pdf(file="MHCIboxplot.pdf",width=5,height=4)                         
p+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.format")
dev.off()

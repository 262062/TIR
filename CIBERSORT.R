
library(ggpubr)
pFilter=0.05
rt=read.table("CIBERSORT.txt",sep="\t",header=T,row.names=1,check.names=F)    
data=rt[rt[,"P-value"]<0.05,]
data=data[,1:(ncol(rt)-3)]
Type=read.table("Immunity.txt",sep="\t",check.names=F,row.names=1,header=F)
Type=Type[row.names(data),]
colnames(Type)=c("cluster","Subtype")
outTab=data.frame()
data=cbind(data,Type)
for(i in colnames(data[,1:(ncol(data)-2)])){
  rt1=data[,c(i,"Subtype")]
  colnames(rt1)=c("expression","Subtype")
  ksTest<-kruskal.test(expression ~ Subtype, data = rt1)
  pValue=ksTest$p.value
  if(pValue<pFilter){
      outTab=rbind(outTab,cbind(rt1,gene=i))
      print(pValue)
  }
}
write.table(outTab,file="data.txt",sep="\t",row.names=F,quote=F)
data=read.table("data.txt",sep="\t",header=T,check.names=F)       
data$Subtype=factor(data$Subtype, levels=c("L","M","H"))
p=ggboxplot(data, x="gene", y="expression", color = "Subtype", orientation = "horizontal",
     ylab="Fraction",
     xlab="",
     palette = c("chartreuse4","blue","red") )
p=p+rotate_x_text(45)
pdf(file="ciberboxplot.pdf",width=12,height=7)                      
p+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.format")
dev.off()

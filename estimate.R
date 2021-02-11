
library(pheatmap)                   
rt=read.table("ssgsea.txt",sep="\t",header=T,row.names=1,check.names=F)   
Type=read.table("Immunity.txt",sep="\t",check.names=F,row.names=1,header=F)
rt=rt[,row.names(Type)]
score=read.table("scores.txt",sep="\t",check.names=F,row.names=1,header=T)
score=score[row.names(Type),]
colnames(Type)=c("cluster","Subtype")
cluster=cbind(Type,score)
cluster=cluster[,-1]
pdf("estimate.pdf",height=5,width=9)
pheatmap(rt, annotation=cluster, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3)
dev.off()


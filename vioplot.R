
library(ggpubr)
Type=read.table("Immunity.txt",sep="\t",check.names=F,row.names=1,header=F)
Type=Type[order(Type[,2]),]
score=read.table("scores.txt",sep="\t",check.names=F,row.names=1,header=T)
score=score[row.names(Type),]
colnames(Type)=c("cluster","Subtype")
cluster=cbind(Type,score)
cluster=cluster[,-1]
cluster$Subtype=factor(cluster$Subtype, levels=c("L","M","H"))
my_comparisons=list(c("L","M"),c("_M","H"),c("H","L"))
pdf(file="ESTIMATEscore.pdf",width=6,height=5)
ggviolin(cluster, x="Subtype", y="ESTIMATEscore", fill ="Subtype", cex.main=4, cex.lab=4, cex.axis=4,
         palette = c("springgreen2", "steelblue2", "indianred1"), 
         add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons, size=4)
dev.off()

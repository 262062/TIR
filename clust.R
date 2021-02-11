


library(sparcl)                                                        
                
data=read.table("ssgsea.txt",sep="\t",header=T,check.names=F,row.names=1)   

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]

hc = hclust(dist(t(data)))

y=cutree(hc,3)              

write.table(y,file="cluster.txt",sep="\t",quote=F,col.names=F)

pdf(file="cluster.pdf",width=50,height=20)
ColorDendrogram(hc, y = y, labels = names(y), branchlength = 0.3,xlab=" ",sub=" ",main = " ")
dev.off()


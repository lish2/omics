g <- as.matrix(read.table("aug_f1.txt",head=T,sep="\t"))
col <- as.numeric(as.factor(g[,1]))+1

d <- as.matrix(read.table("aug_f1.txt",head=T,sep="\t")[,3:10])
rownames(d) <- g[,2]

### red: clinical index
### green: drugs
### blue: host propty

par(mfrow=c(2,2))
dd <- d[,1:2]
a <- dd[,2] < 0.05; o <- order(dd[a,1])
barplot(dd[a,1][o],col=col[a][o],xlab="effect size",horiz=T,las=1,cex.names=.8)

dd <- d[,3:4]
a <- dd[,2] < 0.05; o <- order(dd[a,1])
barplot(dd[a,1][o],col=col[a][o],xlab="effect size",horiz=T,las=1,cex.names=.8)

dd <- d[,5:6]
a <- dd[,2] < 0.05; o <- order(dd[a,1])
barplot(dd[a,1][o],col=col[a][o],xlab="effect size",horiz=T,las=1,cex.names=.8)

dd <- d[,7:8]
a <- dd[,2] < 0.05; o <- order(dd[a,1])
barplot(dd[a,1][o],col=col[a][o],xlab="effect size",horiz=T,las=1,cex.names=.8)

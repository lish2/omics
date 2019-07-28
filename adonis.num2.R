library(fdrtool)
library(ade4)
library(vegan)

d <- as.matrix(read.table("BM_CKD.txt",head=T,sep="\t",row.names=1))
g <- as.matrix(read.table("index_ckd.txt",head=T,sep="\t",row.names=1))
###########################################

g <- as.numeric(g[,21])

out <- matrix(0,ncol=3,nrow=ncol(d))
rownames(out) <- colnames(d)
colnames(out) <- c("corr","p","q")

for(i in 1:ncol(d)){
	out[i,1] <- cor(d[,i],g)
	out[i,2] <- cor.test(d[,i],g)$p.value 
}
out[,3] <- qvalue(out[,2])$q
write.table(out,"f",sep="\t",quote=F,append=F)

###########################################
d2 <- as.matrix(read.table("module_CON.txt",head=T,sep="\t",row.names=1))
data <- rbind(d,d2)
#g <- c(rep("D",218),rep("N",66))
g <- c(rep("D",223),rep("N",69))
adonis(data~g, permutations = 999, method = "bray")

d <- as.matrix(read.table("MGS_ckd.txt",head=T,sep="\t",row.names=1))
g <- as.matrix(read.table("table.txt",head=T,sep="\t",row.names=1))
g <- g[rownames(d),]

out <- matrix(0,ncol=2,nrow=ncol(g))
rownames(out) <- colnames(g)

for(i in 1:ncol(g)){
	a <- adonis(d~as.factor(g[,i]), permutations = 999, method = "bray")
	out[i,1] <- a$aov.tab[1,5]
	out[i,2] <- a$aov.tab[1,6]
}

adonis(d~as.factor(g[,1])*as.factor(g[,15]), permutations = 999, method = "bray")
adonis(d~as.factor(g[,13])*as.factor(g[,15]), permutations = 999, method = "bray")
adonis(d~as.factor(g[,1])*as.factor(g[,13]), permutations = 999, method = "bray")
adonis(d~as.factor(g[,1])*as.factor(g[,13])*as.factor(g[,15]), permutations = 999, method = "bray")
adonis(d~as.factor(g[,13])*as.factor(g[,14]), permutations = 999, method = "bray")
adonis(d~as.factor(g[,13])*as.factor(g[,21]), permutations = 999, method = "bray")

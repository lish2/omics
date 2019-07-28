library(reshape2)
library("fdrtool")
#library(circlize)

num = function(x){
  m = seq(-1,1,0.1)
  for(i in m ){
    x[i<x&x<(i+0.1)]<-i
  }
  return(x)
}
data_p = as.matrix(read.table("corr_mgs_all_cor_spearman_p",sep = "\t",check.names = F,
                              header = T,row.names = 1,stringsAsFactors = F))
data_r = as.matrix(read.table("corr_mgs_all_cor_spearman_r",sep = "\t",check.names = F,
                                  header = T,row.names = 1,stringsAsFactors = F))
for(i in 1:ncol(data_r)){
  for(j in i:col(data_r)){
    data_r[i,j] = 0
  }
}

data_r_he = melt(data_r)
data_p_he = melt(data_p)
data_p_he = data_p_he[which(data_r_he[,3]!=0),]
data_r_he = data_r_he[which(data_r_he[,3]!=0),]

BM = read.table("BM_CKD.txt",check.names = F,sep ="\t",header = T,row.names = 1)
FM= read.table("FM_CKD.txt",check.names = F,sep = "\t",header = T,row.names = 1)
L7 = read.table("mgs_ckd.txt",check.names = F,sep = "\t",header = T,row.names = 1)
#Med = read.table("Med_CKD.txt",check.names = F,sep= "\t",header = T,row.names = 1)

data_r_he$type = "positive"
data_r_he$type[data_r_he[,3]<0] = "negative"
data_r_he$num = num(data_r_he[,3])
#data_r_he[,1] = gsub("[|]p__.*","",data_r_he[,1])
#data_r_he[,2] = gsub("[|]p__.*","",data_r_he[,2])
data_r_he$one = "L7"
data_r_he$two = "L7"

data_r_he[as.character(data_r_he[,1])%in%as.character(colnames(BM)),][,6] = "BM"
data_r_he[as.character(data_r_he[,2])%in%as.character(colnames(BM)),][,7] = "BM"
data_r_he[as.character(data_r_he[,1])%in%as.character(colnames(FM)),][,6] = "FM"
data_r_he[as.character(data_r_he[,2])%in%as.character(colnames(FM)),][,7] = "FM"
data_r_he$all = paste(data_r_he[,6],data_r_he[,7],sep = "_")
data_p_he = data_p_he[data_r_he$one!=data_r_he$two,]
data_r_he = data_r_he[data_r_he$one!=data_r_he$two,]

adf_p = fdrtool(data_p_he$value,statistic = "pvalue")
data_p_he$qvalue = adf_p$qval
data_r_he = data_r_he[data_p_he$qvalue<=0.01,]
data_p_he = data_p_he[data_p_he$qvalue<=0.01,]
range(abs(data_r_he$value))
data_p_he = data_p_he[abs(data_r_he[,3])>=0.35,]
data_r_he = data_r_he[abs(data_r_he[,3])>=0.35,]

a = c()
for(i in 1:length(unique(data_r_he[,2]))){
  b = sum(unique(data_r_he[,2])[i]==data_r_he[,2])
  a = c(a,b)
}

write.table(data_r_he,"edge_network_mgs_CKD",quote = F,sep = "\t",row.names = F)

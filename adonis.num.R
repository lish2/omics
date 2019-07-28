library(vegan)
library(psych)
###########Get rid of the correlation > 0.5#########
corr_self <- function(x,y){
    cutoff = 0.5
    name = x[,1]
    MGS1 = y[,name]
    corr_MGS_FM = corr.test(MGS1)
    a = corr_MGS_FM$r
    for(i in 1:nrow(a)){
        for(j in i:nrow(a)){
            a[j,i] = 0
        }
    }
    for(k in 1:nrow(x)){
        a = a[!a[k,]>cutoff,!a[k,]>cutoff]
        if(nrow(a)==k){
            break
        }
    }
    return(row.names(a))
}
##################################################################

############The selected variables are calculated adonis##########
calculta_R <- function(x,y){
    x_FM <- x[,c(1,2,3)]
    x_BM <- x[,c(1,4,5)]
    x_FM <- x_FM[order(x_FM[,2],decreasing = T),]
    x_BM <- x_BM[order(x_BM[,2],decreasing = T),]
    x_FM <- x_FM[x_FM$P<=0.05,]
    x_BM <- x_BM[x_BM$P<=0.05,]
    y_FM <- y[,corr_self(x_FM,y)]
    y_BM <- y[,corr_self(x_BM,y)]
    y_FM_adonis <- adonis(FM~.,data = y_FM)
    y_BM_adonis <- adonis(BM~.,data = y_BM)
    a = list(FM = y_FM_adonis,BM = y_BM_adonis)
    a
}
################file###############################
FM_name = "FM_CON.txt"
BM_name = "BM_CON.txt"
type1_name = "mgs_con.txt"#host_prop_ckd.txt
type2_name = "module_con.txt"#index_ckd.txt
name1 = "mgs"#host
name2 = "module"#index
all_name = paste(name1,name2,sep = "_")
write_name = "Classified_information_CKD"#clinical_CKD

###############read file#########################
FM = read.table(FM_name,header = T,row.names = 1,sep = "\t",
                check.names = F,stringsAsFactors = F)
BM = read.table(BM_name,header = T,row.names = 1,sep = "\t",
                check.names = F,stringsAsFactors = F)
host_prop = read.table(type1_name,header = T,sep = "\t",row.names = 1,
                 check.names = F,stringsAsFactors = F)
index = read.table(type2_name,header = T,sep = "\t",row.names = 1,
                 check.names = F,stringsAsFactors = F)
index = index[row.names(host_prop),]
#index[,42] = as.double(index[,42])
index <- index[,!apply(index,2,function(x){sum(is.na(x))/nrow(index)>0.5})]
for(i in 1:ncol(index)){
    index[,i][which(is.na(index[,i]))] = mean(index[,i],na.rm = T)
}
host_prop <- host_prop[,!apply(host_prop,2,function(x){sum(is.na(x))/nrow(host_prop)>0.5})]
for(i in 1:ncol(host_prop)){
    host_prop[,i][which(is.na(host_prop[,i]))] = mean(host_prop[,i],na.rm = T)
}###################NA-> mean

#MGS = index
all_two = cbind(index,host_prop)

adonis_index = read.table("module_adonis_CKD.txt",sep = '\t',check.names = F,
                           stringsAsFactors = F,header = T)#index_all_CKD_index.txt
adonis_host = read.table("mgs_adonis_CKD.txt",sep = '\t',check.names = F,
                          stringsAsFactors = F,header = T)#host_prop_all_CKD_index.txt
#adonis_mgs = adonis_index
adonis_mgs = rbind(adonis_index,adonis_host)

index_adonis = calculta_R(adonis_index,index)
index_adonis_FM = index_adonis$FM
index_adonis_BM = index_adonis$BM
index_adonis_FM_adj = RsquareAdj(1-index_adonis_FM[[1]][5][,1][length(index_adonis_FM[[1]][5][,1])-1],
                                 66,length(index_adonis_FM[[1]][5][,1])-2)
index_adonis_BM_adj = RsquareAdj(1-index_adonis_BM[[1]][5][,1][length(index_adonis_BM[[1]][5][,1])-1],
                                 66,length(index_adonis_BM[[1]][5][,1])-2)

host_adonis = calculta_R(adonis_host,host_prop)
host_adonis_FM = host_adonis$FM
host_adonis_BM = host_adonis$BM
host_adonis_FM_adj = RsquareAdj(1-host_adonis_FM[[1]][5][,1][length(host_adonis_FM[[1]][5][,1])-1],
                                 66,length(host_adonis_FM[[1]][5][,1])-2)
host_adonis_BM_adj = RsquareAdj(1-host_adonis_BM[[1]][5][,1][length(host_adonis_BM[[1]][5][,1])-1],
                                 66,length(host_adonis_BM[[1]][5][,1])-2)


all_index_host = calculta_R(adonis_mgs,all_two)
all_FM = all_index_host$FM
all_BM = all_index_host$BM
all_FM_adj = RsquareAdj(1-all_FM[[1]][5][,1][length(all_FM[[1]][5][,1])-1],
                                66,length(all_FM[[1]][5][,1])-2)
all_BM_adj = RsquareAdj(1-all_BM[[1]][5][,1][length(all_BM[[1]][5][,1])-1],
                                66,length(all_BM[[1]][5][,1])-2)

data1 = matrix(rep(0,24),ncol = 4)
colnames(data1) = c("type1","type2","value","number")
data1[,1] = c(name1,name1,name2,name2,all_name,all_name)
data1[,2] = c(rep(c("BM","FM"),3))
data1[,3] = c(host_adonis_BM_adj,host_adonis_FM_adj,index_adonis_BM_adj,index_adonis_FM_adj,
              all_BM_adj,all_FM_adj)
data1[,4] = c(length(host_adonis_BM[[1]][5][,1])-2,length(host_adonis_FM[[1]][5][,1])-2,
              length(index_adonis_BM[[1]][5][,1])-2,length(index_adonis_FM[[1]][5][,1])-2,
              length(all_BM[[1]][5][,1])-2,length(all_FM[[1]][5][,1])-2)

write.table(data1,write_name,quote = F,sep = "\t",row.names = F)

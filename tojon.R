library(jsonlite)
library(RJSONIO)
data1 = read.table("edge_network_mgs_CKD",sep = "\t",header = T,stringsAsFactors = F)
data1$type[data1$type=="positive"] = "#4c91c1"
data1$type[data1$type=="negative"] = "#ff0000"
data1[,1] = paste(data1[,6],data1[,1],sep = ".")
data1[,2] = paste(data1[,7],data1[,2],sep = ".")

a = c(unique(c(data1[,1],data1[,2])))
color_type = c()
json_file1 = vector("list",length(a))
w = 0 
for(i in 1:length(a)){
    if(a[i]%in%data1[,1]){
        b = data1[data1[,1]==a[i],]
        c = b[,2]
        color = b[,4]
    }else{
        b = data1[data1[,2]==a[i],]
        c = b[,1]
        color = b[,4]
    }
    color_type = c(color_type,color)
    json_file1[[i]]$name = a[i]
    json_file1[[i]]$size = sample(10000,1)
    if(length(c)==1){
        json_file1[[i]]$imports =list(c)
    }else{
        json_file1[[i]]$imports =c
    }
    
}

k = paste(color_type,collapse = "\",\"")
#grep(k,"\\","")
write.csv(k,"color.txt",row.names = F)

json_file2 = toJSON(json_file1)
writeLines(json_file2, "finCKD_mgs.json")

cat(toJSON(json_file1))

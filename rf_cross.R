suppressPackageStartupMessages(library(randomForest))
library(tcltk)  
suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser()
parser$add_argument("-n", help = "input number")

args <- parser$parse_args()

infile <- file.path(args$i)
innum <- file.path(args$n)

cat("\nUsing the following setting:\n")
cat("Input classified: ", innum, "\n")

Species = read.table("mgs.profile",row.names = 1,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
FM = read.table("BM.txt",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
#FM = data.frame(t(FM),check.names = F)
#for(i in 1:nrow(FM)){
#  FM[i,][which.max(FM[i,])] = mean(as.double(FM[i,]))
#}
#apply(FM,2,function(x){x[which.max(x)]=mean(x)})
FM = sweep(FM, 2, apply(FM,2,sum), "/")
Species = sweep(Species,2,apply(Species,2,sum),"/")
#FM = FM[-grep("unknown",row.names(FM)),]
#FM = FM[c("Phenol","Indole","Phenylacetaldehyde","p-Cresol","Acetic acid","Butyric acid","Propionic acid"),]
#FM = data.frame(t(FM),check.names = F)
#Species = data.frame(t(Species),check.names = F)
Species = Species[row.names(FM),]
bb = colnames(FM)
a = bb[as.numeric(innum)] 
FM = FM[,a]
#Species = Species[!is.na(FM),]
#FM = FM[!is.na(FM)]
#FM[which(FM==max(FM))]  = max(FM[-which(FM==max(FM))])
names(Species) <- make.names(names(Species))
fit <- randomForest(FM~ ., data=Species, importance=TRUE, proximity=TRUE, ntree=1000)
imp <- fit$importance
impvar <- imp[order(imp[,1],decreasing = TRUE),]
Species = Species[,row.names(impvar)]
RMSE = matrix(rep(0,nrow(Species)*30),nrow = nrow(Species),ncol = 30)
colnames(RMSE) = colnames(Species)[1:30]
row.names(RMSE) = row.names(Species)
set.seed(123)
for(j in 2:30){
  Species1 = Species[,1:j]
  for(i in 1:nrow(Species)){
    traindata = Species1[-i,]
    train_FM = FM[-i]
    testdata = Species1[i,]
    test_FM = FM[-i]
    fit <- randomForest(train_FM~ ., data=traindata, importance=TRUE, proximity=TRUE, ntree=1000)
    preds <- predict(fit, testdata)
    RMSE[i,j] <- as.double(preds)
  }
  print(j)
}
write.table(RMSE,paste(a,"predicted.txt",sep = "_"),quote = F,sep = "\t")

result = matrix(rep(0,150),nrow = 30,ncol = 5)
colnames(result) = c("RMSE","DRMSE",'COR','Q2',"name")
result[,5] = colnames(Species)[1:30]
for(m in 2:30){
  rmse = sqrt(mean((as.double(FM)-as.double(RMSE[,m]))^2))
  drmse = rmse/(max(as.double(FM))-min(as.double(RMSE[,m])))
  Q2 = 1-sum((FM-as.double(RMSE[,m]))**2)/sum((FM-mean(as.double(RMSE[,m])))**2)
  result[m,3] = cor(as.double(FM),as.double(RMSE[,m]))
  result[m,1] = rmse 
  result[m,2] = drmse
  result[m,4] = Q2
}
write.table(result,paste(a,"result.txt",sep = "_"),quote = F,sep = "\t")

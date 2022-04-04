library(singleCellNet)
source("/fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_TrainingDataSelection/Scripts/functions.R")


args = commandArgs(trailingOnly=TRUE)
print(args)
folder <- args[1]
tag <- args[2]
testfolder <- args[3]
out <- args[4]
print("Starting....")
input <- paste(sep="/", folder, tag)
output <- paste(sep="/", out,tag)
print(output)
set.seed(100) #can be any random seed number

data <- readData(paste(sep="/", input, "data_train.txt"))
meta <- read.csv(paste(sep="/", input, "meta_train.txt"))


test <- readData(paste(sep="/", testfolder, "data_test.txt"))
meta_test <- read.csv(paste(sep="/", testfolder, "meta_test.txt")) 
print(head(meta_test))

ncells <- as.numeric(stringr::str_split(tag,pattern="_", simplify = T)[2])
stList = splitCommon(sampTab=meta, ncells=ncells, dLevel="class_")
meta = stList[[1]]
data = data[,meta$id]

print("-------------------")
class_info<-scn_train(stTrain = meta, expTrain = as.matrix(data), nTopGenes = 10, nRand = 100, nTrees = 1000, 
                      nTopGenePairs = 25, dLevel = "class_", colName_samp = "id")

print("-------------------")
nqRand = 50
predictions <- as.data.frame(t(scn_predict(class_info[['cnProc']], test, nrand=nqRand)))

predicted <- unlist(colnames(predictions)[apply(predictions,1,which.max)])


predictions$id = rownames(predictions)
predictions$predicted = predicted


results <- merge(predictions, meta_test, by="id", all = T)
print(head(results))
results$prediction.match <- results$predicted == results$class_
write.table(results, paste0(output,".txt") ,sep="\t", quote=F)
source("/fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_TrainingDataSelection/Scripts/functions.R")


getExperiment <- function(datafile, metafile=NULL){
   data <- readData(datafile)
   experiment <- SingleCellExperiment::SingleCellExperiment(list(counts=data))
   if(!is.null(metafile)){
     meta <- read.csv(metafile)
     SingleCellExperiment::colLabels(experiment) <- meta$class_
     experiment <- experiment[,!is.na(experiment$label)]   
   }
  experiment <- scater::logNormCounts(experiment)
  return(experiment)
}
args = commandArgs(trailingOnly=TRUE)
folder <- args[1]
tag <- args[2]
testfile <- args[3]
out <- args[4]

input <- paste(sep="/", folder, tag)
output <- paste(sep="/", out,tag)

train <- getExperiment(paste(sep="/", input, "data_train.txt"), paste(sep="/", input, "meta_train.txt"))
test <- getExperiment(paste(sep="/", testfile, "data_test.txt"))

predictions<- as.data.frame(SingleR::SingleR(test=test, ref=train, labels=train$label, de.method="wilcox"))
predictions$id <- rownames(predictions)

meta_test <- read.csv(paste(sep="/", testfile, "meta_test.txt"))
results <- merge(meta_test, predictions, by="id")
colnames(results)[colnames(results)== "labels"] <- "predicted"
results$prediction.match <- results$predicted == results$class_

write.table(results, paste0(output,".txt") ,sep="\t", quote=F)

if (!require("Seurat")) install.packages("Seurat")
library("Seurat")

readData <- function(data){
  df <-  as.data.frame(data.table::fread(data,sep = ","))
  rownames(df) <- df$V1
  df[,1] <- NULL
  return(df)
}

prepareData <- function(data, meta, filtering=T, hvg=500, getHVG=T){
  counts <- as.data.frame(readData(data))

  if(filtering){
    counts <- counts[,!(colSums(counts)  * 0.2 < 
                          colSums(counts[startsWith(rownames(counts),"MT-"),]))]
    counts <- counts[,!(colSums(counts != 0) < 200)]
    counts <- counts[,!(colSums(counts) < 1000)]
    counts <- counts[!(rowSums(counts != 0) <10),]
  }
  metadata <- read.csv(meta)
  metadata <- metadata[metadata$id %in% colnames(counts),]
  counts <- counts[, colnames(counts) %in% metadata$id]
  rownames(metadata)<- metadata$id
  obj <-  CreateSeuratObject(counts = counts, meta.data = metadata) 
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = hvg, verbose = FALSE)
  return(obj)
}

predictCelltype_seurat <- function(dataTrain, dataTest, metaTrain, metaTest,
                                   output, path){
  print(dataTrain)
  train <- prepareData(dataTrain, metaTrain, filtering=F)
  print(dataTest)
  test <- prepareData(dataTest, metaTest, filtering=F)
  
  print("Scale data")
  train <- ScaleData(train, verbose = FALSE)
   
  print("Get anchors")
  target.anchors <- FindTransferAnchors(reference = train, query = test,
                                        dims = 1:30, verbose = F)
  print("Transfer data")
  predictions <- TransferData(anchorset = target.anchors,
                              refdata = train$class_, dims = 1:30, verbose = F)
  predictions$predicted.id    <- test$id
  predictions$prediction.match <- predictions$predicted.id == predictions$class_
  write.table(predictions, paste0( path,"_predictions.txt") ,sep="\t", quote=F)
  #target.copy <- test
  #target.copy <- AddMetaData(target.copy, metadata = predictions)
  #target.copy$prediction.match <- target.copy$predicted.id == target.copy$class_
  
  #print("Write data")
    
  #data <- data.frame(target.copy$id, target.copy$class_, target.copy$predicted.id,
  #                       target.copy$prediction.match)
  #colnames(data) <- c("id", "class_", "predicted", "prediction.match")
  #write.table(data, paste0(path,".txt") ,sep="\t", quote=F)
  print("Done")
}


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

input <- args[1]
test <- args[2]

predictCelltype_seurat(paste(sep="/", input, "data_train.txt"),
                       paste(sep="/", test, "data_test.txt"),
                       paste(sep="/", input, "meta_train.txt"),
                       paste(sep="/", test, "meta_test.txt"),
                       basename(input), args[3])


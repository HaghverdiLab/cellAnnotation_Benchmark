#library(Seurat)
preprocessData <- function(data, tfs=NULL, hvg=NULL, filtering=F, genes=NULL){
  if(filtering){
    data <- data[,!(colSums(data)  * 0.2 < 
                      colSums(data[startsWith(rownames(data),"MT-"),]))]
    data <- data[,!(colSums(data != 0) < 200)]
    data <- data[,!(colSums(data) < 1000)]
    data <- data[!(rowSums(data != 0) <10),]
  }
  if(!is.null(genes)) data <- data[rownames(data) %in% genes,]
  data_seurat <- Seurat::CreateSeuratObject(counts=data)
  data_seurat <- Seurat::NormalizeData(data_seurat, verbose = FALSE)
  if(!is.null(hvg)){
    data_seurat <- Seurat::FindVariableFeatures(data_seurat, selection.method = "vst",
                                        nfeatures = hvg, verbose = FALSE)
    data <- Seurat::GetAssayData(object = data_seurat, slot = "data")
    if(is.null(tfs)) data <- data[rownames(data) %in%
                                    Seurat::VariableFeatures(data_seurat),]
    else data <- data[rownames(data) %in% Seurat::VariableFeatures(data_seurat)|
                        rownames(data) %in% tfs,]
  }else data <- Seurat::GetAssayData(object = data_seurat, slot = "data")
  
  data <- data.matrix(data)
  return(data)
}

readData <- function(data){
  df <-  as.data.frame(data.table::fread(data,sep = ","))
  rownames(df) <- df$V1
  df[,1] <- NULL
  return(df)
}
getMeta <- function(metafile, presentTypes=NULL){
  meta <- read.csv(metafile)
  meta <- meta[order(meta$class_),]
  if(!is.null(presentTypes)){
    meta <- rbind(meta[meta$class_ %in% presentTypes,],
                  meta[!(meta$class_ %in% presentTypes),] )
  }
  return(meta)
}
getFileName <- function(folder, type, set){
  name <- paste(folder, paste(type, set, sep="_"), sep="/")
  return(paste(name, "txt", sep="."))
}

getDataFile <- function(data, output, meta, set){
  meta <- meta[meta$id %in% colnames(data),]
  data <- data[, meta$id]
  write.table(data, getFileName(output, "data", set), sep=",", quote=F,
              row.names = F, col.names = T)
  write.table(meta, getFileName(output, "classes", set), sep=",", quote=F,
              row.names = F, col.names = T)
  return(data)
}

args = commandArgs(trailingOnly=TRUE)

folder <- args[1]
tag <- args[2]
test <- args[3]
out <- args[4]

input <- paste(sep="/", folder, tag)
output <- paste(sep="/", out,tag)
if(!dir.exists(output))dir.create(output)

dataTrain <- readData(paste(sep="/", input, "data_train.txt"))
dataTest <-readData(paste(sep="/", test, "data_test.txt"))
meta_test <- getMeta(paste(sep="/", test, "meta_test.txt"))
meta_train <-  getMeta(paste(sep="/", input, "meta_train.txt"), unique(meta_test$class_))
dataTrain <- preprocessData(dataTrain, hvg=500)
dataTest <- preprocessData(dataTest, hvg=NULL, genes = rownames(dataTrain))
trainingData <- getDataFile(dataTrain,output, meta_train, "train")
testData <- getDataFile(dataTest, output, meta_test,"test")
write.table(rownames(trainingData), paste(output, "feature_names.txt", sep="/"),
              sep=",", quote=F, row.names = F, col.names = F)


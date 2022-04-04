getMeta <- function(metafile, presentTypes=NULL){
  meta <- read.csv(metafile)
  meta <- meta[order(meta$class_),]
  if(!is.null(presentTypes)){
    meta <- rbind(meta[meta$class_ %in% presentTypes,],
                  meta[!(meta$class_ %in% presentTypes),] )
  }
  return(meta)
}

getDataFile <- function(data, output, meta, set){
  meta <- meta[meta$id %in% rownames(data),]
  data <- data[meta$id,]
  if(!file.exists(output))dir.create(output) 
  write.table(data, getFileName(output, "data", set), sep=",", quote=F,
              row.names = T, col.names = T)
  write.table(meta, getFileName(output, "meta", set), sep=",", quote=F,
              row.names = F, col.names = T)
}

getFileName <- function(folder, type, set){
  name <- paste(folder, paste(type, set, sep="_"), sep="/")
  return(paste(name, "txt", sep="."))
}

readData <- function(data){
  df <-  as.data.frame(data.table::fread(data,sep = ",", verbose = F))
  rownames(df) <- df$V1
  df[,1] <- NULL
  return(df)
}

print("Starting....")
args = commandArgs(trailingOnly=TRUE)
tag <-args[1]
print(tag)
test <- args[3]
#test <- paste(sep="/", args[3], tag)
train <- paste(sep="/", args[2], tag)

output <- paste(sep="/", args[4], tag) 

print("Read..")

data_test <- t(readData(paste(sep="/",test, "data_test.txt")))
data_train <- t(readData(paste(sep="/",train, "data_train.txt")))


print("Get meta...")
meta_test <- getMeta(paste(sep="/",test, "meta_test.txt"))
print(unique(meta_test$class_))
meta_train <- getMeta(paste(sep="/",train, "meta_train.txt"), unique(meta_test$class_))

print("Get data...")
trainingData <- getDataFile(data_test,output, meta_test,"test")


#if(stringr::str_detect(tag, "full")){
#data_train <- t(read.table(paste(sep="/",train, "data_test.txt"), sep=",", header=T))
#meta_train <- getMeta(paste(sep="/",train, "meta_test.txt"))
trainingData <- getDataFile(data_train,output, meta_train,"train")
#}

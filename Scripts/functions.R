######## General Functions ######## 
 
#' Read the expression data matrix (fast)
#' Uses fread to increase the reading spead of the expression data
#' @param data file containing the expression matrix´
readData <- function(data){
  df <-  as.data.frame(data.table::fread(data,sep = ",", verbose = F))
  rownames(df) <- df$V1
  df[,1] <- NULL
  return(df)
}

## Highly variable gene selection ##
pearson_residuals <- function(counts, theta){
    counts_sum1 = rowSums(counts)
    counts_sum0 = colSums(counts)
    counts_sum  = sum(counts)

    #get residuals
    mu = (counts_sum1  %*% t(counts_sum0)) / counts_sum
    z = (counts - mu) / sqrt(mu + mu**2/theta)

    #clip to sqrt(n)
    n = ncol(counts)
    z[z >  sqrt(n)] = sqrt(n)
    z[z < -sqrt(n)] = -sqrt(n)
    return(as.matrix(z))
}

select_hvg <- function(data, hvgs){
  residuals = pearson_residuals(data,200)
  residual_var = matrixStats::rowVars(residuals)
  names(residual_var) <- rownames(residuals)
  residual_var <- names(residual_var[order(residual_var, decreasing = T)])
  data = data[head(residual_var, hvgs), ]
  return(data)
}


## Preprocessing  the data ##
prepareData <- function(data, meta, project, hvgs=200, features=NULL){
    counts <- as.data.frame(readData(data))
    metadata <- read.csv(meta)
    metadata <- metadata[metadata$id %in% colnames(counts),]
    counts <- counts[, colnames(counts) %in% metadata$id]
    rownames(metadata)<- metadata$id
    
    counts <- preprocessing(counts, features=features, hvgs=hvgs)
    obj <-  CreateSeuratObject(counts = counts, meta.data = metadata) 

  return(obj)
}

preprocessing <- function(data, features=NULL, hvgs=200){
    if(!is.null(features)){
      data <- data[rownames(data) %in% features, ]
    } else data <- select_hvg(data, hvgs)
    data <- log2(data + 0.001)
    data <- as.matrix(data)
    data <- data/colSums(data)[col(data)]
    return(data)
}

getFileName <- function(folder, type, set){
  name <- paste(folder, paste(type, set, sep="_"), sep="/")
  return(paste(name, "csv", sep="."))
}

getFiles <- function(datafile, metafile=NULL,set=NULL, output=NULL, features=NULL, hvgs=500){
    print(paste(set, output))
    data <- readData(datafile)
    metadata <- read.csv(metafile)
    metadata <- metadata[metadata$id %in% colnames(data),]
    metadata <- metadata[order(metadata$class_),]
    data <- data[, colnames(data) %in% metadata$id]
    rownames(metadata)<- metadata$id
    if(!is.null(features))  data <- data[rownames(data) %in% features, ]
    # We skip the preprocessing step since ItClust has internal processing
    #data <- preprocessing(data, features=features, hvgs=hvgs)
    data <- data[order(rownames(data)),]
    
    metadata <- metadata[metadata$id %in% colnames(data),]
    data <- data[, metadata$id]

    write.table(t(data), getFileName(output, "data", set), sep=",", quote=F,
              row.names = F, col.names = T)
    
    write.table(colnames(data), getFileName(output, "cells", set), sep=",", quote=F,
                          row.names = F, col.names = F)

    write.table(metadata, getFileName(output, "meta", set), sep=",", quote=F,
              row.names = F, col.names = T)
    return(data)
}

getExperiment <- function(datafile, metafile=NULL, features=NULL, hvgs=200){
    data <- readData(datafile)
    metadata <- read.csv(metafile)
    metadata <- metadata[metadata$id %in% colnames(data),]
    data <- data[,colnames(data) %in% metadata$id]
    rownames(metadata)<- metadata$id
    data <- preprocessing(data, features=features, hvgs=hvgs)
    data <- data[, order(colnames(data))]

    experiment <- SingleCellExperiment::SingleCellExperiment(list(logcounts=data))
    metadata <- metadata[metadata$id %in% colnames(data),]
    SingleCellExperiment::colLabels(experiment) <- metadata$class_[order(metadata$id)]

  return(experiment)
}



getAccuracy <- function(data, col){
    pred <- data$class_ == data[, col]
    pred[pred == TRUE] <- 1
    pred[is.na(pred)] <- 0
    return( sum(pred) / nrow(data))
}
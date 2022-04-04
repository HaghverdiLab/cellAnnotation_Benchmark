library(docstring)
######## General Functions ######## 
 
#' Read the expression data matrix (fast)
#' 
#' Uses fread to increase the reading spead of the expression data
#' @param data file containing the expression matrix´
readData <- function(data){
  df <-  as.data.frame(data.table::fread(data,sep = ",", verbose = F))
  rownames(df) <- df$V1
  df[,1] <- NULL
  return(df)
}


getDE <- function(c1, c2,output){
    if(c1 == c2) return(list(c1,c2,NULL))
    df <- FindMarkers(pbmc, ident.1 = c1, ident.2 = c2, min.pct = 0.5)
    DE <- rownames(df)
    c1 <- unlist(stringr::str_split(c1, " "))[1]
    c2 <- unlist(stringr::str_split(c2, " "))[1]
    print(paste(c1,c2))
    write.table(DE,
        paste(output, stringr::str_replace_all(c1, " ", "_"), stringr::str_replace(c2, " ", "_"), sep="_"), quote=F, row.names=F)
    return(list(c1,c2, DE))
}


getMeans <- function(d, datasets){
   genes <- unlist(lapply(datasets, function(set) length(d$x[d$x %in% set$V1])))
   return(mean(genes))                      
}

add_column_DE <- function(folder, pattern, method, df, name){
    folders <- list.files(folder, full.names = T, pattern = pattern)
    dd <- folders[stringr::str_detect(folders,method)]
    datasets <- lapply(dd, function(folder) 
                       as.data.frame(read.csv(paste(sep="/", folder, "feature_names.txt"),
                                              header=F), header=F))
    df[,name] <- (unlist(lapply(data, function(d) getMeans(d, datasets)))) 
    return(df)                           
}                          
     
                                
summarize_data <- function(data, label=NULL, steps=NULL){
    if(!is.null(label)) data <- data[data$Version == label,]
    if(!is.null(label)) data <- data[data$CellsPerCelltype %in% steps,]
    
    data_melt <- reshape2::melt(data, value.name = "Accuracy",
                                id.vars = c("Method", "Version", "CellsPerCelltype", "SetNr"))
    
    colnames(data_melt)[colnames(data_melt) == "variable"] <- "CellType"
    data_summary <- data_melt %>%
                    group_by(Method, Version, CellsPerCelltype, CellType) %>%
                    summarise(mean = mean(Accuracy),
                              sd= sd(Accuracy),
                              q25= quantile(Accuracy, probs=0.25),
                              q75= quantile(Accuracy, probs=0.75),
                              median = median(Accuracy))
    
    return(as.data.frame(data_summary))   
}

get_percentage_predictions <- function(predictions,version, method, size, nrSets=20){
    name <- paste(sep="_", version, method)
    set <-  paste(sep="_", version, method, size)
    predictions[, name] <- rowSums(predictions[, stringr::str_detect(colnames(predictions), set)])
    predictions[, name] <- (predictions[, name] * 100) / nrSets
    predictions[is.na(predictions[,name]), name ] <- 0 
    return(predictions)
}

prepare_umap <- function(file_train, file_test, meta_data,split){
 data_train <- readData(file_train)
 data_test  <- readData(file_test)
 data <- cbind(data_train, data_test)
 
 pbmc <- Seurat::CreateSeuratObject(data, meta.data = meta_data)
 print(pbmc)
 if(!is.null(split)){
    pbmc.list <- SplitObject(pbmc, split.by = split)
     for (i in 1:length(pbmc.list)) {
    pbmc.list[[i]] <- NormalizeData(pbmc.list[[i]], verbose = FALSE)
    pbmc.list[[i]] <- FindVariableFeatures(pbmc.list[[i]], selection.method = "vst",
                                           nfeatures = 2000, verbose = FALSE)
    }
    pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, dims = 1:10, verbose = F)
    pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, dims = 1:10, verbose = FALSE, k.weight=10)
    DefaultAssay(pbmc.integrated) <- "integrated"
 } else{
    pbmc.integrated <- NormalizeData(pbmc, verbose = FALSE)
    pbmc.integrated  <- FindVariableFeatures(pbmc.integrated, selection.method = "vst",
                                           nfeatures = 2000, verbose = FALSE)
 }
 pbmc.integrated  <- Seurat::ScaleData(pbmc.integrated, verbose = F)
 pbmc.integrated  <- Seurat::RunPCA(pbmc.integrated, features = Seurat::VariableFeatures(object = pbmc.integrated), verbose = F)
 pbmc.integrated  <- Seurat::RunUMAP(pbmc.integrated, dims = 1:10, verbose = F) 
 return(pbmc.integrated)
}

getData <- function(folder, attention_file, seed, attentionNames, expressionNames){
    folder <- paste(sep="_", folder, seed)
    attention_file <- paste0(attention_file, "_",seed, "/attentionmap.txt")
    
    train <- read.csv(paste(sep="/", folder, "data_train.txt"))
    features <- read.csv(paste(sep="/", folder, "feature_names.txt"), header=F)
    rownames(train)<- features$V1
    
    meta <- read.csv(paste(sep="/", folder, "classes_train.txt"))
    types <- unique(meta$class_)
    
    classes <- lapply(types, function(type) train[,colnames(train) %in% meta$id[meta$class_ == type]])
    means <-  lapply(classes, function(class) rowMeans(class))
    names(means) <- types
                        
    train_summary <- do.call(cbind, means)
    colnames(train_summary)<- expressionNames
 
    test <- read.csv(paste(sep="/", folder, "data_test.txt"))
    rownames(test)<- features$V1
    
    meta <- read.csv(paste(sep="/", folder, "classes_test.txt"))
    types <- unique(meta$class_)               
    
    classes <- lapply(types, function(type) test[,colnames(test) %in% meta$id[meta$class_ == type]])
    means <-  lapply(classes, function(class) rowMeans(class))
    names(means) <- types
    
    test_summary <- do.call(cbind, means)
    colnames(test_summary)<- paste0("Test", expressionNames)
                       
                     
    attention <- read.csv(attention_file, sep="\t")
    attention <- attention[,2:ncol(attention)]
    colnames(attention) <- attentionNames
    rownames(attention) <- features$V1

    data <- cbind(train_summary,attention)
    data <- cbind(test_summary,data)
                     
    corr <- cor(data)
    corr_lower <- get_lower_tri(corr)
    cor_melt <- reshape2::melt(corr_lower, na.rm = TRUE)

    cor_melt <- tidyr::separate(cor_melt,Var1,into = c("Type1","Celltype1"),sep = " ",remove = TRUE,extra = "merge")
    cor_melt <- tidyr::separate(cor_melt,Var2,into = c("Type2","Celltype2"),sep = " ",remove = TRUE,extra = "merge")
    cor_melt$Type1<- factor(cor_melt$Type1, levels = c("Expression", "Attention", "TestExpression"))
    cor_melt$Type2<- factor(cor_melt$Type2, levels = c("Expression", "Attention", "TestExpression"))
    cor_melt$Seed <- seed
    return(cor_melt)
}

                 
# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }


### Method to read the results from the result folder for one of the methods
# Input:
#  resultfolder - folder including the results from the run_scripts.R script
#  id           - Identifyer for the reference dataset used
#  method       - Method used to generate the results
#
# Output:  Dataframe containing a summary of the results in the result folder
get_results_method <- function(resultfolder, id, method){
    print(paste("Start", method, "..."))
    files <- list.files(resultfolder, pattern=id, full.names = T)
    if (method == "Seurat") files <- files[stringr::str_detect(files,
                                                               "predict",
                                                               negate=T)]
    data <- lapply(files, function(file) getVector(method, file))
    summary <- data %>% reduce(full_join, by = "id")
    return(summary)
}

get_results <- function(folder, name, methods, query){
    datasets <- lapply(methods,
                   function(method) get_results_method(paste(sep = "/",
                                                             folder,
                                                             method),
                                                       name, method) )
                    
    data <- datasets %>% reduce(full_join, by = "id")                   
    data <- merge(data,query) 
    data <- adjust_names(data) # make sure cell types have the same naming conventions
    rownames(data) <- data$idprint
    return(data)    
}

### Get the vector containing the result summary for one result file
# Input:
#  method - method used for the predictions
#  folder - file / folder containing the results
#
# Output: Dataframe containing for each celll the ground truth,
#         the predicted label and the cell id
getVector <- function(method, folder){
    if(method %in% c("Seurat", "SingleCellNet", "CellID", "SingleR")){
        if(stringr::str_detect(folder, "predict")) stop("Wrong file")
        name <- paste(sep="_", getName(folder,1),
                      method , getName(folder,2),
                      stringr::str_replace(getName(folder,3), ".txt", ""))
  df <- read.csv(folder, sep="\t")
  df <- df[,c("id", "predicted")]
  colnames(df) <- c("id", name)
    
} else if (method == "ItClust"){
  name <- paste(sep="_",  getName(folder,1), "ItClust",getName(folder,2),
                getName(folder,3))
  df <- read.csv(paste(folder, "results.txt", sep="/"), sep="\t")
  df <- df[,c("class_", "predicted_celltype", "cell_id")]
  colnames(df) <- c("class", "predicted", "id")
  df[,name] <- df$predicted 
  df <- df[, c("id",  name)]
}
    rownames(df) <- df$id
    return(df)
} 
                   
### Function to get part of the identifyer of a dataset
# Input:
#  folder   - name of the file/folder containing the results
#  position - position of the information needed
#
# Output: Infomtion about the parameters for the result file
getName <- function(folder, position){
  name <- unlist(stringr::str_split(folder,"/"))
  name <- name[length(name)]
  name <- unlist(stringr::str_split(name,"_"))[position]
  return(name)
}
                   
### Adjust the cell labels in a subset of columns
# Input:
#  data - dataset containing the predictions
#  name - identifyer for the columns of interest
# Output: Dataframe with the translated columns
adjust_names <- function(data, name="ItClust"){
    cols <- colnames(data)[stringr::str_detect(colnames(data), name)] 
    x<- do.call(cbind,lapply(cols, function(col) translate(data[,col])))
    colnames(x) <- cols  
    id <- data[,!(colnames(data) %in% cols)]

    data <- cbind(x,id)
    return(as.data.frame(data))
}

### Helper functions with the translations that need to be made
# Input: Column that needs to be translated
#
# If working with a different dataset this is the function that
# That needs to be adjusted
translate <- function(col){
    col[col == "B"] <- "B cell"
    col[col == "Cytotoxic T"] <- "Cytotoxic T cell"
    col[col == "CD4+ T"] <- "CD4+ T cell"
    col[col == "Dendritic"] <- "Dendritic cell"
    col[col == "Natural killer"] <- "Natural killer cell"
    col[col == "Plasmacytoid dendritic"] <- "Plasmacytoid dendritic cell"
    return(col)
} 
                            
### Turn the results into long fromat
# Input:
#  data      - result dataframe
#  celltypes - list of cell types in the data (factor)
#  methods   - list of methods in the data (factor)
#  sizes     - number of cells in the full reference data for each celltype
#  cols      - columns next to the prediciton column included in the data
#
# Output: Results in long format
make_long <- function(data, celltypes, methods, sizes=NULL,
                      cols= c("id","nGene", "nUMI", "percent.mito",
                              "Cluster", "class_", "Experiment",
                              "Method")){
    x <- reshape2::melt(data,  id.vars = cols,
                        value.name = "Prediction")
    x[c('Reference','Approach',
        "Size", "Set")] <- stringr::str_split_fixed(x$variable, '_', 4)
    
    x$Match <- x$Prediction == x$class_ 
    if(!(is.null(sizes)))x$refSize <- sizes[ match(x$class_, celltypes ) ]
    x$Approach <- factor(x$Approach, levels=methods)
    x$class <- factor(x$class_, levels=celltypes)
    x <- x[,c("id", "Prediction", "Reference", "Approach", "Size", "Set",
              "class", "refSize", "Match")]
    x$Size <- as.numeric(x$Size)
    x <- x[!is.na(x$Match),]
    return(x)
}  
                             
### Get the dataset with the result summary
# Input:
#  data - results in long format
#  type - celltype of interest
#  method - method of interest
#  set -  set of interest
get_measures <- function(data, type, method,size=NULL, set=NULL){
    data <- data[data$method == method,]
    if(!is.null(size)) data <- data[data$size == size,]
    if(!is.null(set)) data <- data[data$set == set,]
    
    tp <- length(data$predicted[data$predicted == type & data$class_ == type])
    fp <- length(data$predicted[data$predicted == type & data$class_ != type])
    fn <- length(data$predicted[data$predicted != type & data$class_ == type])
    tn <- length(data$predicted[data$predicted != type & data$class_ != type])
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    f1 <- 2*(precision * recall) / (precision + recall)
    accuracy <- (tp) / length(data$predicted[data$class_ == type])
    return(data.frame("class"=type,"method"=method,"set"=set,
                      "precision"=precision,"recall"=recall,"f1"=f1,
                      "accuracy"=accuracy))
}
                             
                          
get_summary_bootstrap <- function(folder, reference, outfolder, query, name,
                                  celltypes, methods, seqs, full){
    data <- get_results(folder, name, methods, query)
    data <- make_long(data, celltypes, methods)
    data <- weighted_bootstrap(read.csv(ref), data)
    data_set <- do.call(rbind, lapply(celltypes, function(type)
                do.call(rbind, lapply(methods,   function(method) 
                do.call(rbind, lapply(seqs, function(set) 
                        get_measures(mono, type, method, set)))))))
    data_set <-  merge(data_set, full, by=c("class", "method"))
    data_id <- summarize_by_id(data) 
                                  
    write.table(data_set,
                paste(sep="/", outfolder,paste0("summary_", name, ".csv" )),
                col.names=T, row.names=T, quote=T, sep=",")
    return(data_id)
}
                                      
summarize_by_id <- function(data){
    data$match <- data$class_ == data$predicted
    data$match[data$match == TRUE] <- 1
    summary <- data %>% 
           dplyr::group_by(id, class_, method) %>% 
           dplyr::summarize(score = sum(match) / n()) 
    return(summary)
}                                      
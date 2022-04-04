library(docstring)

######## Data preparation - Construction of Training and Validation sets ######## 

#' Match expression matrix and meta data
#' 
#' Reduces the expression matrix to cells that are also included in the meta data
#' @param data expression matrix with cells in the columns
#' @param meta meta data for the expression matrix 
#' @param out_data file in which to save the resulting expression matrix
#' @param out_meta file in which to save the resulting meta data
get_subset <- function(data, meta, out_data, out_meta){
  
  data <- data[, colnames(data) %in% meta$id]
  print(paste(ncol(data), nrow(data)))
  
  write.table(data, out_data, row.names = T,
              col.names = T, quote = F, sep=",")
  print(nrow(meta))
  meta <- meta[meta$id %in% colnames(data),]
  write.table(meta, out_meta, row.names = F,
              col.names = T, quote = F, sep=",")
  print(nrow(meta))
}


#' Creation of a set of random sets
#' 
#' Creates a set of random sets based on different numbers of maximum cells per
#' celltype and seed for each set.
#' @param input input folder containing file data_train.txt and meta_train.txt
#' @param output folder in which to write the results
#' @param n_random number of sets oer size. Default=1
#' @param name identifier of the datasets, Default=current date
#' @param steps sizes for which to construct sets. Default= c(100,500,1000,2000,3000)
get_random_sets <- function(input, output, n_random=1, name=NULL, steps=NULL){
 if(is.null(name)) name <- stringr::str_replace_all(paste(Sys.Date()), "-","")
 if(is.null(steps)) steps <- append(seq(1000, 3000, 1000), c(500,100))
 if(!dir.exists(output)) dir.create(output) 
    
 meta <- read.csv(paste(sep="/", input, "meta_train.txt"))
 data <- readData(paste(sep="/", input, "data_train.txt")) 

 ids <- Reduce(intersect,list(meta$id, colnames(data)))
 meta <- meta[meta$id  %in% ids,]    
   
 sets <- seq(1,n_random,1)

 x <- lapply(steps,
             function(step) lapply(sets,
                                   function(set) write_random(data, meta, step,
                                                             set,name, output)))
}
                                          
write_random <- function(data, meta, step, set, name, output){
    out <- paste(sep="/", output, paste(sep="_",name, step, set))
    if(!dir.exists(out)) dir.create(out)
    
    meta_sub <- do.call(rbind,
                        lapply(unique(meta$class_),
                               function(type) get_random(meta, step, type,set)))
    
    write.table(meta_sub, paste(sep="/", out, "meta_train.txt"), sep=",",
              row.names = F, col.names = T, quote = F)
  
    data_sub <- data[,colnames(data) %in% meta_sub$id]
    write.table(data_sub, paste(sep="/", out, "data_train.txt"), sep=",",
              row.names = T, col.names = T, quote = F)  
}           

get_random <- function(data, nr ,type, seed){
  df <- data[data$class_ == type,]
  if(nrow(df)<=nr)return(df)  
  set.seed(seed)
  df <- df[sample(nrow(df), nr),]
  return(df)
}                               
                               
get_validation_set <- function(full, meta, train, output){
  data_train <- read.csv(train)
  meta <- meta[!(meta$id %in% data_train$id),]
  data <- full[,(colnames(full) %in% meta$id)]
  if(!dir.exists(output))dir.create(output)
  write.table(meta, paste(output, "meta_test.txt", sep="/"), col.names = T,
              row.names = T, sep=",", quote = F)
  write.table(as.matrix(data), paste(output, "data_test.txt", sep="/"),
              col.names = T, row.names = T, sep=",", quote = F)
}

addX <- function(name){
  if (startsWith(name, "4")) return(paste0("X",name))
  else return(name)
}
updateID <- function(idVector){
  idVector <- gsub("-", ".", idVector)
  idVector <- unlist(lapply(idVector, function(n) addX(n)))
  return(idVector)
}
getMetaFormat <- function(meta){
  meta <- meta[,c("celltype", "index")]
  colnames(meta) <- c("class_", "id")
  meta$id <- updateID(meta$id)
  return(meta)
}

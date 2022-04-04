#!/usr/bin/env Rscript

visualizeItClust <- function(folder){
  clustering <- read.csv(paste(folder, "results/clustering_results.csv", sep="/"))
  ids <- read.csv(paste(folder, "results/celltype_assignment.txt", sep="/"), header=F)
  
 # if(stringr::str_detect(folder, "PBMC")){
  #  truth <- read.csv("/fast/AG_Haghverdi/Carla_Moelbert/itClust/PBMC_full/meta_test.txt", sep=",")
  #}else truth <- read.csv("/fast/AG_Haghverdi/Carla_Moelbert/itClust/CrossSpecies_full/meta_test.txt", sep=",")
  

  truth <- read.csv(paste(folder, "meta_test.txt", sep="/"), sep=",")
  ids <- strsplit(as.character(ids$V1), " be ")
  ids <- as.data.frame(do.call(rbind,ids))
  ids$V1 <- gsub("Cluster ","",ids$V1)
  ids$V2 <- gsub(" cell","",ids$V2)
  
  x <- as.data.frame(do.call(rbind, strsplit(ids$V1, " ")))
  data <- data.frame(ids$V2, x$V1)
  colnames(data) <- c("class_", "cluster")
  clustering$predicted_celltype <- data$class_[ match(clustering$cluster,
                                                      data$cluster ) ]
  
  clustering$cell_id <- gsub("-target","", clustering$cell_id)
  x <- merge(clustering, truth, by.x="cell_id", by.y = "id", )
  x$class_ <- gsub(" cell","", x$class_)
  name <- paste(sep="/", folder, "results.txt")
  write.table(x, name, col.names = T, row.names = T, quote = F, sep="\t")
}

args = commandArgs(trailingOnly=TRUE)

folder <- args[1]
visualizeItClust(folder)

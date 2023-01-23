source("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Notebooks/functions.r")
                                                                   
args = commandArgs(trailingOnly=TRUE)
folder <- args[1]
name <- args[2]
sizes <- unlist(stringr::str_split(pattern=",", args[3]))
meta <- read.csv(args[4]) 
  
x <- expand.grid(list(name, sizes, seq(1,20,1)))

tags_all <- list.files(paste(sep="/", folder, "ItClust"), pattern=name)

tags <- apply(x, 1, function(row) paste(collapse ="_", unlist(row), sep=""))
tags <- stringr::str_replace_all(tags, " ", "")

tags <- tags[tags %in% tags_all]
             

              
              
seurat <- do.call(rbind,lapply(tags, function(tag) readConfidence_Seurat(paste(sep="/", folder, "Seurat"), tag)))
                             
scn <- do.call(rbind,lapply(tags, function(tag) readConfidence_SCN(paste(sep="/", folder, "SingleCellNet"), tag)))
                        
singleR <- do.call(rbind,lapply(tags, function(tag) readConfidence_SingleR(paste(sep="/", folder, "SingleR"), tag)))
                              
itclust <- do.call(rbind,lapply(tags, function(tag) readConfidence_Itclust(paste(sep="/", folder, "ItClust"), tag,
                                                                meta)))
                                
data <- rbind(itclust, seurat, scn, singleR)
print(head(data))                               
print                                
write.csv(data, args[5], col.names = T, row.names=F, sep=",", quote=F)
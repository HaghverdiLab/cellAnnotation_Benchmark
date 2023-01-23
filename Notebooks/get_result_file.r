setwd(getSrcDirectory()[1])

source("functions.R")
source("visulizations.r")
args = commandArgs(trailingOnly=TRUE)
folder <- args[1]
metafile <- args[2]
name <- args[3]                        
out <- args[4]
out_long <- args[5]
 
    
                             


                             
                             
write.table(data, out, sep=",", col.names=T, row.names=T, quote=F, append=F)  
 


long <- transform_PBMC_results(out, celltypes, methods, sizes)  

                             
all <- long[long$Set <= 20 & long$Size %in% c(38,100,250, 500,1000,1500,2000,3000),]
measures <- do.call(rbind, lapply(unique(data$class),
                                  function(type) do.call(rbind,lapply(unique(data$Approach),
                                  function(method) do.call(rbind, lapply(unique(data$Size),
                                  function(size) do.call(rbind, lapply(unique(data$Set),
                                  function (set) get_measures(all, type, "PBMC10x",
                                                              method, size,set)))))))))
                                                              
write.csv(measures, "../Results/Files/values_all.csv")  
                                                                       
                                                                       
write.table(data, out_long, sep=",", col.names=T, row.names=T, quote=F, append=F)
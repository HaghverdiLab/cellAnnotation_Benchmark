setwd("/fast/AG_Haghverdi/Carla_Moelbert/cellAnnotation_Benchmark/Notebooks/")
library(ggplot2)
library(dplyr)
library(purrr)
library(Seurat)
library(viridis)
library("RColorBrewer")
source("../Scripts/functions.R")
source("../Scripts/prepare_results_functions.r")

##########################################################################

celltypes = c("Cytotoxic T cell", "CD4+ T cell", "CD14+ monocyte", 
              "B cell", "Megakaryocyte", "Natural killer cell",
              "CD16+ monocyte", "Dendritic cell",
              "Plasmacytoid dendritic cell")

methods <- c("Seurat", "SingleR","CellID", "SingleCellNet", "ItClust")  
sizes <- c(3090, 2418, 1373, 1022, 703, 623, 273, 126, 38)
names(sizes) <- celltypes

query <- read.csv("../../Celltype_annotation/Data/Fulldata/PBMC_Query/meta.csv")
folder <- "../../Celltype_annotation/Data/Predictions/"
name <- "PBMC10x"
number_of_sets <- 20 
dataset_sizes <- c(38,100,250, 500,1000,1500,2000,3000)
##########################################################################

# Write a daframe with the results from each method
data <- get_results(folder, name, methods, query) 

long <- make_long(data, celltypes, methods, sizes) 

print("---------------------------")
full <- do.call(rbind, lapply(unique(long$class), function(type)
        do.call(rbind,lapply(unique(long$Approach), function(method)
                get_measures(long, type, method, size=3090, set=0)))))
print(head(full))

print("---------------------------")                                                                        
all <- long[long$Set <= number_of_sets & long$Size %in% dataset_sizes,]     
all <- do.call(rbind, lapply(unique(long$class), function(type)
       do.call(rbind,lapply(unique(long$Approach), function(method) 
       do.call(rbind, lapply(unique(long$Size), function(size) 
       do.call(rbind, lapply(unique(long$Set), function (set) 
               get_measures(data=all, type=type, method, size,set)))))))))
                                                                       
print(head(all))                                                             
write.csv(all, "../Results/Files/values_all.csv")  
write.csv(full, "../Results/Files/values_full.csv")   
rm(all)

print("---------------------------") 
full <- full[, c("class", "method", "precision", "recall", "f1", "accuracy")]
colnames(full) <-  c("class", "method", "full_precision", "full_recall", "full_f1",
                     "full_accuracy")   
                            
mono <- get_summary_bootstrap("../../Celltype_annotation/Data/Predictions_abundancebased//",
                    "../Data/Fulldata/PBMC10x_Reference/meta.csv",
                    "../Results/Files/", "Mono", celltypes, methods,
                    seq(1,20,1), full)
print("---------------------------")
mosaic <- get_summary_bootstrap("../../Celltype_annotation/Data/Predictions_mosaic/",
                    "../Data/Fulldata/PBMC10x_Reference/meta.csv",
                    "../Results/Files/", "Mosaic", celltypes, methods,
                      seq(1,20,1), full)
                            
print("---------------------------")
                            
full <- data[, stringr::str_detect(colnames(data), "3090")]
full$tech <- data$Method
full <- reshape2::melt(full,id=c("class_", "id", "tech"),
                       value.name = "predicted")
full$score <-  full$class_ == full$predicted
full$score[full$score == TRUE] <- 1
full[c('reference', 'method',
       "size", "set")] <- stringr::str_split_fixed(full$variable, '_', 4)
full <- full[, c("id", "class_", "method", "score", "predicted", "tech")]
mosaic <- tidyr::pivot_wider(mosaic_id, names_from = c(method), values_from = score,
                             names_pref="mosaic_")
mono <- tidyr::pivot_wider(mono_id, names_from = c(method), values_from = score,
                             names_pref="mono_")

full1 <- tidyr::pivot_wider(full[, colnames(full) != "predicted"], names_from = c(method),
                            values_from = score,
                             names_pref="full_")

full2 <- tidyr::pivot_wider(full[, colnames(full) != "score"], names_from = c(method),
                            values_from = predicted,
                             names_pref="fullPred_")


data <- Reduce(function(x, y) merge(x, y, all=TRUE), list(mosaic, mono, full1, full2))
write.table(data, "../Results/Files/umap_data.csv", sep=",")                            
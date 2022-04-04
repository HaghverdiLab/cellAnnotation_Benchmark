BiocManager::install("CelliD")
library("CelliD")
source("/fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_TrainingDataSelection/Scripts/functions.R")
getSeuratObject <- function(datafile, metafile, project){
    data <- readData(datafile)
    meta <- read.csv(metafile)
    rownames(meta) <- meta$id
    obj <- CreateSeuratObject(counts = data, project = project, min.cells = 5, meta.data = meta)
    return(obj)
}

args = commandArgs(trailingOnly=TRUE)
folder <- args[1]
tag <- args[2]
testfile <- args[3]
out <- args[4]

input <- paste(sep="/", folder, tag)
output <- paste(sep="/", out,tag)

print("Prepare training data")
print(input)
train <- getSeuratObject(paste(sep="/", input, "data_train.txt"),paste(sep="/", input, "meta_train.txt"), "train")

# Library-size normalization, log-transformation, and centering and scaling of gene expression values
train <- NormalizeData(train, verbose = F)
train <- ScaleData(train, features = rownames(train), verbose = F)

train <- RunMCA(train)
train <- RunPCA(train, features = rownames(train), verbose = F)
train <- RunUMAP(train, dims = 1:30, verbose = F)

train_cell_gs <- GetCellGeneSet(train, dims = 1:50, n.features = 200)
train_group_gs <- GetGroupGeneSet(train, dims = 1:50, n.features = 200, group.by = "class_")



print("Prepare test data")
test <- getSeuratObject(paste(sep="/", testfile, "data_test.txt"),paste(sep="/", testfile, "meta_test.txt"), "test")

test <- NormalizeData(test, verbose = F)
test <- FindVariableFeatures(test, verbose = F)
test <- ScaleData(test, verbose = F)
test <- RunMCA(test, nmcs = 50)

print("Cell type prediction")
predicted_cell_gs <- RunCellHGT(test, pathways = train_cell_gs, dims = 1:50)
predicted_cell_gs_match <- rownames(predicted_cell_gs)[apply(predicted_cell_gs, 2, which.max)]
predicted_cell_gs_prediction <- train$class_[predicted_cell_gs_match]
predicted_cell_gs_prediction_signif <- ifelse(apply(predicted_cell_gs, 2, max)>2, yes =predicted_cell_gs_prediction, "unassigned")

print(head(predicted_cell_gs_prediction_signif))
test$predicted <- predicted_cell_gs_prediction_signif


print("Save results")
test$prediction.match <- test$predicted == test$class_

data <- data.frame(test$id, test$class_, test$predicted,
                         test$prediction.match)
colnames(data) <- c("id","class_", "predicted", "prediction.match")

write.table(data, paste0(output,".txt") ,sep="\t", quote=F)
print("Done")
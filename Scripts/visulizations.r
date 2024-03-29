prepare_umap <- function(file,  meta_data, normalization_method="L1", min_cells=3, hvgs=500, split=NULL){
    set.seed(0)
    data <- readData(file)
  
    #data <- data[, colnames(data) %in% rownames(meta_data)]
    meta_data <- meta_data[rownames(meta_data)%in% colnames(data),]
    data <- data[!(rowSums(data != 0) < min_cells),]

    meta_data <- meta_data[meta_data$id %in% colnames(data),]
    rownames(meta_data) <- meta_data$id
    print(hvgs)
    data <- select_hvg(data, hvgs)
    if(!is.null(split)){
      sets <- lapply(unique(meta_data[, split]), 
                     function(s) data[, colnames(data) %in%   meta_data$id[meta_data[,split] == s]])
      data <- batchelor::fastMNN(sets)
      print("....")
      data <-  SummarizedExperiment::assay(data, "reconstructed")
    } else{

        # Cell wise (Result: Sum of each row = 1)
        if(normalization_method == "L1"){
            data <- apply(data, 2, normalize_l1)
        }else if(normalization_method == "L2") data <- apply(data, 2, normalize_l2)

    }

    print("Scaling...")
    # Gene wise (Variance of each column is 1, mean = 0) (Colum should be Gene)
    data <- t(scale(t(data), center = TRUE, scale = TRUE))
    
    print("create seurat object....")   
    pbmc <- Seurat::CreateSeuratObject(as.matrix(data), meta.data = meta_data)

    print("Set assay....")                 
    pbmc <- Seurat::SetAssayData(pbmc, assay = "RNA", slot = "scale.data", new.data = as.matrix(data))
    
    print("PCA....")
    pbmc <- Seurat::RunPCA(pbmc, features = rownames(data),verbose = F)
    print("Get UMAP...")
    pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10, verbose = F)
    return(pbmc)
}

                     getVectors <- function(method, data){
    summary <- data[data$method == method,]  %>% 
           dplyr::group_by(id) %>% 
           dplyr::summarize("accuracy_confidence" = mean(accuracy_confidence),
                            "accuracy_prediction" = mean(accuracy_prediction),
                            "accuracy_cf" = mean(accuracy_confidenceFiltered),
                            "accuracy_pf" = mean(accuracy_predictionFiltered ),
                           "accuracy_full" = mean(accuracy_full))
    colnames(summary) <- c("id", paste(sep="_", "accuracy_confidence", method) , paste(sep="_","accuracy_prediction", method),
                           paste(sep="_","accuracy_cf", method) , paste(sep="_","accuracy_pf", method), paste(sep="_","accuracy_full", method))
    return(summary)
}


translate <- function(col){
  
    col[col == "B"] <- "B cell"
    col[col == "Cytotoxic T"] <- "Cytotoxic T cell"
    col[col == "CD4+ T"] <- "CD4+ T cell"
    col[col == "Dendritic"] <- "Dendritic cell"
    col[col == "Natural killer"] <- "Natural killer cell"
    col[col == "Plasmacytoid dendritic"] <- "Plasmacytoid dendritic cell"
    return(col)
}     
                     

transform_PBMC_results <- function(data_file,  celltypes, methods, sizes){
    data <- read.csv(data_file)

    x <- reshape2::melt(data,  id.vars =c("id"),value.name = "Prediction")
    
    x[c('Reference', 'Approach', "Size", "Set")] <- stringr::str_split_fixed(x$variable, '_', 4)
    x$Match <- x$Prediction == x$class_ 

    x$refSize <- sizes[ match(x$class_, celltypes ) ]
    x$Approach <- factor(x$Approach, levels=methods)
    x$class <- factor(x$class_, levels=celltypes)
    x <- x[,c("id", "Prediction", "Reference", "Approach", "Size", "Set", "class", "refSize", "Match")]
    x$Size <- as.numeric(x$Size)

    return(x)
}

calculate_accuracy <- function(data){
    summary <- data %>% 
           dplyr::group_by(Size, Approach, Reference, class, Set, refSize) %>% 
           dplyr::summarize(accuracy= mean(Match)) 

    summary <- summary[order(summary$refSize), ]
    return(summary)
}

plot_features <- function(data,groups, split=NULL,nrow=1, ncol=1, cells=NULL, return_plotlist=FALSE,
                          legend_ncol=3, title=""){

    umaps <- lapply(groups, function(method) Seurat::FeaturePlot(data, reduction="umap", 
                                                                 split.by=split,
                                                                 features=method) +
                    labs(title=title, color= "Accuracy")+ 
                    theme(axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks=element_blank(),
                          plot.title=element_text(size=8),
                          legend.text=element_text(size=8), legend.title=element_text(size=8),
                          legend.key.size = unit(0, 'lines'),
                         legend.justification = "center")+
                    guides(color = guide_legend(ncol = legend_ncol, override.aes = list(size=3)))+
                    viridis::scale_color_viridis(guide = "colourbar")
                   ) 
                   
    names(umaps) <- groups
    if(return_plotlist)return(umaps)
    #plot <- ggpubr::ggarrange(plotlist=umaps, common.legend = F, nrow=nrow, ncol=ncol)
    return(umaps)
}                    
                     

                        
plot_umap <- function(data,groups, split=NULL,color="Spectral",
                      nrow=1, ncol=1, cells=NULL, return_plotlist=FALSE, legend_ncol=3, title=""){

    umaps <- lapply(groups, function(method) Seurat::DimPlot(data,reduction="umap", split.by=split, 
                                                     group.by=method, 
                                                     cols=color[names(color) %in% data@meta.data[,method]],
                                                     cells = cells,label=F)+
                   labs(title=title)+ 
                    theme(axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks=element_blank(),
                          legend.text=element_text(size=8), legend.position="right", 
                          plot.title=element_text(size=8),
                          legend.key.size = unit(0, 'lines'))+
                    guides(color = guide_legend(ncol = legend_ncol, override.aes = list(size=3))) 
                   ) 
                   
    names(umaps) <- groups
    if(return_plotlist)return(umaps)
    #plot <- ggpubr::ggarrange(plotlist=umaps, common.legend = F, nrow=nrow, ncol=ncol)
    return(umaps)
}

get_data_lineplot <- function(data, maxsize=3000){
    summary <- data %>% 
           dplyr::group_by(size, method, reference, class) %>% 
           dplyr::summarize(mean_accuracy = mean(accuracy),
                            p25_accuracy = quantile(accuracy, probs = c(0.25)),
                            p75_accuracy = quantile(accuracy, probs = c(0.75)),
                            mean_precision = mean(precision),
                            p25_precision = quantile(precision, probs = c(0.25)),
                            p75_precision = quantile(precision, probs = c(0.75))) 
    
    return(summary)
}
                        
                        
normalize_l1<- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}                        
                        
                        
                        
getVectors <- function(method, data){
    summary <- data[data$method == method,]  %>% 
           dplyr::group_by(id) %>% 
           dplyr::summarize("accuracy_confidence" = mean(accuracy_confidence),
                            "accuracy_prediction" = mean(accuracy_prediction),
                            "accuracy_cf" = mean(accuracy_confidenceFiltered),
                            "accuracy_pf" = mean(accuracy_predictionFiltered ),
                           "accuracy_full" = mean(accuracy_full))
    colnames(summary) <- c("id", paste(sep="_", "accuracy_confidence", method) , paste(sep="_","accuracy_prediction", method),
                           paste(sep="_","accuracy_cf", method) , paste(sep="_","accuracy_pf", method), paste(sep="_","accuracy_full", method))
    return(summary)
}

getUmapData <- function(file, metafile, methods){
    data_bootstrap <- read.csv(file)
    data_bootstrap <- data_bootstrap  %>% 
           dplyr::group_by(id, method, class_) %>% 
           dplyr::summarize(accuracy_full= mean(match_full),
                            accuracy_confidence = mean(match_confidence),
                            accuracy_prediction = mean(match_binary),
                            accuracy_confidenceFiltered = mean(match_cf),
                            accuracy_predictionFiltered = mean(match_filtered)) 
    meta =read.csv(metafile)
    bootstrap_sets <- lapply(methods[methods != "CellID"], function(method) getVectors(method, data_bootstrap))
    bootstrap <- Reduce(function(x, y) merge(x, y, by="id"),bootstrap_sets)
    umapdata <- merge(meta, bootstrap, by="id",all.x=TRUE )
    rownames(umapdata)<- umapdata$id
    return(umapdata)                                             
}
#####################################################################################################################

plot_confidence_scores <- function(data){
   plot <- ggplot(data, aes(as.factor(size), score, color=match)) + 
    geom_boxplot(alpha=0.3, width=0.5, position = position_dodge(0.5),  outlier.size=0.1)+
    facet_grid(cols=vars(class), rows=vars(method),  labeller = label_wrap_gen(width=15))+
    theme_bw()+
    theme(legend.position = "bottom", legend.text=element_text(size = 8),
          legend.title=element_text(size = 8), legend.key.size = unit(0, 'lines'),
          axis.text=element_text(size = 8), axis.title=element_text(size = 8),
          strip.text = element_text(size = 8)) +
    xlab("Cells per cell type")+ylab("Confidence Scores")+labs(color="Predictions")+ylim(0,1)+
    scale_y_continuous(breaks=c(0,0.5,1))
    return(plot)
}
                        
get_violin_plot <- function(data, colors, celltypes, methods, title, value, color){
    data$class <- factor(data$class, levels=celltypes)
    data$method <- factor(data$method, levels=methods)
    
    if(value =="Accuracy"){
      plot <- ggplot(data, aes(class, accuracy, color=class, fill=class)) +
              geom_boxplot(alpha=0.3)+
              geom_point(aes(y=full_accuracy), color="black",
                         size= 0.5, show.legend = F)
    } else {
      plot <- ggplot(data, aes(class, precision, color=class, fill=class)) + 
              geom_boxplot(alpha=0.3)+
              geom_point(aes(y=full_precision), color="black", size= 0.5,
                         show.legend = F)
    }
    plot <- plot + 
            facet_wrap(facets=vars(method),labeller = label_wrap_gen(width=12),
                       nrow = 1) + labs(color="", fill="", title=title)
    plot <- addFormatting(plot, value, "Cell type", "bottom", "none", color)
    return(plot)
}

addFormatting <- function(plot, ylab="", xlab="", legend="bottom", xtext="none",
                          celltypes=NULL, ylim = c(0,1), size= 8, keysize = 0.5){
    if(!is.null(celltypes)) plot <- plot + scale_fill_manual(values=celltypes) +
                                           scale_color_manual(values=celltypes)
    plot <- plot + 
            theme_bw() + ylab(ylab) + xlab(xlab) + 
            theme(axis.text=element_text(size=size), axis.title=element_text(size=size),
                 legend.position = legend, legend.key.size = unit(keysize, 'lines'),
                 legend.text=element_text(size = size),
                 strip.text = element_text(size = size), title=element_text(size = size)) +
            ylim(ylim) 
    
    if(xtext=="none"){
        plot <- plot + theme( axis.text.x = element_blank(),
                             axis.ticks.x = element_blank())
    } else if (xtext == "angle") {
        plot <- plot + theme(axis.text.x=element_text(angle=25, hjust = 1))
    }

    return(plot)
}
                        
get_lineplot <- function(data, value){
    
    if(value == "Accuracy"){
        plot = ggplot(data, aes(size, mean_accuracy, group=method)) +
           geom_ribbon(aes(ymin = p25_accuracy, ymax = p75_accuracy),
                       alpha = 0.3)+
           geom_hline(aes(yintercept = accuracy), color="red",
                      linetype="dashed")
    } else {
        plot = ggplot(data, aes(size, mean_precision, group=method)) +
               geom_ribbon(aes(ymin = p25_precision, ymax = p75_precision),
                       alpha = 0.3)+
           geom_hline(aes(yintercept = precision), color="red",
                      linetype="dashed")
    }
    
    plot <- plot + geom_point(size=0.5) + geom_line()+ 
            geom_vline(aes(xintercept=refSize), linetype="dotted",
                       color="black")+
            facet_grid(rows = vars(class), cols =vars(method), 
                       labeller = label_wrap_gen(width=10)) + 
            scale_y_continuous(breaks=c(0,0.5,1))
    
    plot <- addFormatting(plot, value, "Cells per cell type", xtext="")
    
    return(plot)
}

                        
get_plot_comparison <-  function(data, value, colors){
    if(value == "Accuracy"){
        plot <- ggplot(data, aes(class, accuracy, color=class, size=set, #shape=set, 
                                  group=set))
    } else plot <-  ggplot(data, aes(class, precision,color=class, size=set, #shape=set,
                                      group=set))
    
    plot <- plot + geom_point(alpha=0.5) + #size=2
            labs(color="Celltype", size= "Gene set used") + 
            facet_grid(cols = vars(method)) +
    scale_size_manual(values = c("1000 HVGs" = 3, "200 HVGs"=1.5))
     plot <- addFormatting(plot, ylab=value, xlab="Cell types", legend="left", xtext="none", celltypes= colors) 

}

                        
get_data_curated <- function(file, path, method){
    data <- read.csv(paste(path,file, sep="/"), sep="\t")
    data$tag <- stringr::str_replace(file, ".txt", "")
    data$method <- method
    data <- data[, c("id", "predicted", "tag", "method")]
    data[c('reference', 'size', "set")] <- stringr::str_split_fixed(data$tag, '_', 3)
    return(data)
}

get_measures <- function(data, type, ref, method, set){
    data <- data[data$ref == ref & data$method == method & data$set == set,] 
    tp <- length(data$predicted[data$predicted == type & data$class_ == type])
    fp <- length(data$predicted[data$predicted == type & data$class_ != type])
    fn <- length(data$predicted[data$predicted != type & data$class_ == type])
    tn <- length(data$predicted[data$predicted != type & data$class_ != type])
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    f1 <- 2*(precision * recall) / (precision + recall)
    accuracy <- (tp) / length(data$predicted[data$class_ == type])
    return(data.frame("class"=type,"reference"=ref,"method"=method,"set"=set,
                      "precision"=precision,"recall"=recall,"f1"=f1, "accuracy"=accuracy))
}
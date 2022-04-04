# Functions to visualize the results
library(Seurat)
library(ggplot2)

#' Visualize the accuracy distribution
#' 
#' Lineplot showing how the accuracy changes over different numbers of cells 
#' celltype for different methods
#' @param data data frame with q25, q75 and median value for each method for each size
#' @param celltypes vector with celltypes of interest
#' @param ncols number of columns in the plot
#' @param pos legend position. Default=bottom
plot_accuracy <- function(data, celltypes=NULL,
                          methods = c("DataDriven", "PriorKnowledge",
                                      "MLP","ItClust", "Seurat" ),
                          ncols=NULL, pos="bottom", maxSizes=NULL){
    
    data <- as.data.frame(data)
    if(!is.null(celltypes))data <- data[data$CellType %in% celltypes,]
    if(is.null(ncols)) ncols <- sqrt(unique(data$CellType))
    
    plot <- ggplot(data, aes(x=CellsPerCelltype, y=median,
                             color=factor(Method, levels = methods),
                             fill=factor(Method, levels = methods)))+
    geom_line() + geom_point() +
    geom_ribbon(aes(ymin=q25, ymax=q75, color=NULL), alpha=0.2)+
    labs(color="Method", fill="Method")+ xlab("Max. cells per celltype") + 
    ylab("Accuracy") + ylim(0,1) + theme_minimal() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=8),
          plot.title=element_text(size=6, hjust = 0.5),
          legend.title=element_text(size=6), legend.text=element_text(size=6),
          legend.position=pos,legend.key.size = unit(0.5,"line"))+
     geom_vline(aes(xintercept = maxSize), maxSizes, linetype="dashed")
    if(length(celltypes)>1) plot <- plot + facet_wrap(~CellType,ncol=ncols)
    return(plot)
}


plot_features <- function(data, method=NULL, cells, sets=NULL, names=NULL,legend="none"){
    if(is.null(sets)) stop("At least on set needed")
    len <- length(sets)
    if(!(is.null(names)) &len != length(names)) stop("Names must be the same length as the sets")
    if(len == 1){
        plot <- FeaturePlot(data, reduction="umap", features=paste0(sets[1], method),
                           cols=viridis(101), cells = cells) +  labs(color="Accuracy")+
                           theme(axis.text=element_text(size=6), axis.title=element_text(size=6),plot.title = element_text(size=6),
                                 legend.text=element_text(size=6), legend.title = element_text(size=6)) + ggtitle(names[1])
        return(plot)
    }else{
        plots <- lapply(seq(1,len,1), function(n) FeaturePlot(data, reduction="umap", features=paste0(sets[n], method),
                           cols=viridis(101), cells = cells)+  labs(color="Accuracy")+
                           theme(plot.title = element_text(size=6),  axis.text=element_text(size=6), axis.title=element_text(size=6),
                                 legend.text=element_text(size=6), legend.title = element_text(size=6)) + ggtitle(names[n]))
        return(plots)
    }
     
    
}

plot_accuracy_violin <- function(data, label, steps, features, ncol=NULL){
    data <- data[data$Version == label,]
    data <- data[data$CellsPerCelltype %in% steps,]
    data_melt <- reshape2::melt(data,id.vars=c("Method", "Version", "CellsPerCelltype", "SetNr"),
                                value.name = "Accuracy")
    colnames(data_melt)[colnames(data_melt) == "variable"] <- "CellType"
    if(is.null(ncol)) ncol=sqrt(length(unique(data_melt$CellType)))

    plot <-ggplot(data_melt, aes(x=factor(Method, levels = features), y=Accuracy,
                               fill=factor(Method, levels = features),
                               color=factor(Method, levels = features))) + 
        geom_violin(alpha=0.3) + geom_boxplot(width=0.1, alpha=0.5) +
        facet_wrap(facets = ~ CellType, ncol = ncol) +
        xlab(NULL) + ylab("Accuracy") + theme_bw()+ ylim(0,1) + 
        theme( axis.text.y=element_text(size=6), axis.title=element_text(size=8),
              axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.title=element_text(size=10),
              legend.position = "bottom", legend.key.size = unit(0.5,"line")) + labs(fill="Method", color="Method")
    return(plot)

}

plot_umap <- function(data,groups, split=NULL,color="Spectral", nrow=1, ncol=1, cells=NULL, return_plotlist=FALSE, legend_ncol=3){
    
    umaps <- lapply(groups, function(method) DimPlot(data, reduction="umap", split.by=split,
                                                     group.by=method, cols=color, cells = cells)+ 
                    theme(axis.text=element_text(size=6), axis.title=element_text(size=6),
                          plot.title=element_blank(),legend.text=element_text(size=6), legend.position="right",
                          legend.key.size = unit(0, 'lines'))+ guides(color = guide_legend(ncol = legend_ncol, override.aes = list(size=3)))
                   )
    names(umaps) <- groups
    if(return_plotlist)return(umaps)                
    #plot <- ggpubr::ggarrange(plotlist=umaps, common.legend = F, nrow=nrow, ncol=ncol)
    return(umaps)
}


plot_correlation <- function(data, legend="right", names=NULL){
    corr_plot <- ggplot(data = data_mean, aes(x=Celltype1, y=Celltype2, fill = mean))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 20, vjust = 1, size = 6, hjust = 1),
       axis.text.y = element_text(size = 6, hjust = 1),
       legend.text = element_text(size = 6, hjust = 1), legend.title = element_text(size = 6, hjust = 1),
      legend.position=legend, strip.text.x=element_text(size=8), 
      axis.title = element_blank())+
 coord_fixed() + facet_wrap(~interaction(Type1, Type2), labeller = as_labeller(names)) #
return(corr_plot)
} 

plot_properties <- function(data, yaxis){
    plot <- ggplot(d, aes(x=Method, y = data[, yaxis], fill = Method, group =Method)) + 
            geom_violin(alpha=0.3) + geom_boxplot(width=0.1, alpha=0.5)+
            theme_minimal()+ylab(yaxis)+
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                  axis.text.y = element_text(size = 6, hjust = 1),
                  legend.text = element_text(size = 6, hjust = 1),
                  legend.title = element_text(size = 6, hjust = 1),
                  legend.position="right", strip.text.x=element_text(size=8), 
                  axis.title = element_blank())
    return(plot)
}
  
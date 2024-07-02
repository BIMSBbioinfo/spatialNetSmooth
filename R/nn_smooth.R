#' perform nn smoothing
#'
#' @param seu    Seurat object
#' @param genes csv-file indicating which genes to use (optional)
#' @param assay   which assay to use (optional)
#' @param graph which neighbour-graph to use, "nn" or "snn"
#' @param a alpha for smoothing (optional)
#'
#' @return vector with smoothed scores
#' @export
#' @importFrom netSmooth netSmooth


nn_smooth <- function(seu, genes = "Datasets - Ikarus - Gene_lists.csv", assay = "Spatial", a = 0.8, graph = "nn"){
  if(inherits(seu, "Seurat")){
    if (!("gsea_rat_norm" %in% colnames(seu@meta.data))) {
      seu <- gseaCalc(seu, genes, assay)
    }
    graphs <- paste(assay, "_", graph, sep="")
    neighbours <- as.matrix(seu@graphs[[graphs]])
    gsea_score <- seu@meta.data$gsea_rat_norm
    gsea_score <- matrix(gsea_score,ncol = 1)
    rownames(gsea_score) <- rownames(neighbours)
    smoothed <- netSmooth(gsea_score, neighbours, alpha = a)
    return(smoothed)
  }else if(inherits(seu, "VoltRon")){
    if (!("gsea_rat_norm" %in% colnames(Metadata(seu)))) {
      seu <- gseaCalc(seu, genes, assay)
    }
    graph <- as.matrix(as_adjacency_matrix(seu@graph$SNN))
    gsea_score <- Metadata(seu)$gsea_rat_norm
    gsea_score <- as.matrix(gsea_score,ncol = 1)
    rownames(gsea_score) <- rownames(Metadata(seu))
    smoothed <- netSmooth(gsea_score, graph, alpha = a)
    return(smoothed) 
      
  }else{
    stop("Input must be either of type Seurat or VoltRon")
  }
}  

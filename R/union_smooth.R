#' union smoothing
#'
#' @param seu Seurat-Object or VoltRon object
#' @param genes csv-file indicating which genes to use (optional)
#' @param assay which assay to use (optional)
#' @param a alpha for smoothing (optional)
#' @param graph which neighbour-graph to use, "nn" or "snn"
#' @param k distance of cell to neighbors
#'
#' @return vector with smoothed scores
#' @export
#'
#' @examples union_smooth(seu, a = 0.4)
#' @importFrom netSmooth netSmooth
#' @import Seurat
#' @import dplyr
union_smooth <- function(seu, genes = "Datasets - Ikarus - Gene_lists.csv", assay = "Spatial", a = 0.8, graph = "nn", k = 7){
  if(inherits(seu, "Seurat")){
    if (!("gsea_rat_norm" %in% colnames(seu@meta.data))) {
      seu <- gseaCalc(seu, genes, assay)
    }
    graphs <- paste(assay, "_", graph, sep="")
    neighbours <- as.matrix(seu@graphs[[graphs]])
    gsea_score <- seu@meta.data$gsea_rat_norm
    gsea_score <- matrix(gsea_score,ncol = 1)
    rownames(gsea_score) <- rownames(neighbours)
    adj_mat <- adj_matrix(seu, 7)
    union <- (adj_mat | neighbours)*1
    gsea_smoothed <- netSmooth(gsea_score, union, alpha = a)
    return(gsea_smoothed)

  }else if(inherits(seu, "VoltRon")){
    if (!("gsea_rat_norm" %in% colnames(Metadata(seu)))) {
      seu <- gseaCalc(seu, genes, assay)
    }
    adj <- as.matrix(as_adjacency_matrix(seu@graph$radius))
    graph <- as.matrix(as_adjacency_matrix(seu@graph$SNN))
    gsea_score <- Metadata(seu)$gsea_rat_norm
    gsea_score <- as.matrix(gsea_score,ncol = 1)
    rownames(gsea_score) <- rownames(Metadata(seu))
    union <- (adj | graph)*1
    smoothed <- netSmooth(gsea_score, union, alpha = a)
    return(smoothed)
  } 
  else{
    stop("Input must be either of type Seurat or VoltRon")
  }
}

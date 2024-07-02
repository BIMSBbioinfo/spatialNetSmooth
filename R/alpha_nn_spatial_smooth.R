#' alpha_nn_spatial_smooth
#'
#' @param seu a Seurat-Object
#' @param genes csv-file indicating which genes to use (optional)
#' @param assay which assay to use (optional)
#' @param a alpha for smoothing (optional)
#' @param alpha for linear-combination alpha*spatial + (1-alpha)*expression
#' @param graph which neighbour-graph to use, "nn" or "snn"
#'
#' @return vector with smoothed scores
#' @export
#'
#' @examples alpha_nn_spatial_smooth(se, assay= "RNA")
#' @importFrom netSmooth netSmooth
#' @import Seurat
#' @import dplyr

alpha_nn_spatial_smooth <- function(seu, genes = "Datasets - Ikarus - Gene_lists.csv", assay = "Spatial", alpha= 0.4, a = 0.8, graph = "nn"){
  if(inherits(seu, "Seurat")){
    if (!("gsea_rat_norm" %in% colnames(seu@meta.data))) {
      seu <- gseaCalc(seu, genes, assay)
    }
    graphs <- paste(assay, "_", graph, sep="")
    neighbours <- as.matrix(seu@graphs[[graphs]])
    gsea_score <- seu@meta.data$gsea_rat_norm
    gsea_score <- matrix(gsea_score,ncol = 1)
    rownames(gsea_score) <- rownames(neighbours)
    adj_mat <- adj_matrix(seu)
    linear_comb <- alpha*adj_mat + (1-alpha)*neighbours
    gsea_smoothed <- netSmooth(gsea_score, linear_comb, alpha = a)
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
    linear_comb <- alpha*adj + (1-alpha)*graph
    smoothed <- netSmooth(gsea_score, linear_comb, alpha = a)
    return(smoothed)
  }else{
    stop("Input must be either of type Seurat or VoltRon")
  }
}
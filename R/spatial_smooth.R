#' perform spatial smoothing
#' @param seu    Seurat object or VoltRon object
#' @param genes csv-file indicating which genes to use (optional)
#' @param assay   which assay to use (optional)
#' @param a alpha for smoothing (optional)
#' @param k distance of cell to neighbors
#' @return vector with smoothed scores
#' @export
#' @importFrom netSmooth netSmooth
#' @import Seurat
#' @import dplyr

spatial_smooth <- function(seu, genes = "Datasets - Ikarus - Gene_lists.csv", assay = "Spatial", a = 0.8, k = 7){
  if(inherits(seu, "Seurat")){
    if (!("gsea_rat_norm" %in% colnames(seu@meta.data))) {
      seu <- gseaCalc(seu, genes, assay)
        }
    adj_mat <- adj_matrix(seu, k)
    gsea_score <- seu@meta.data$gsea_rat_norm
    gsea_score <- matrix(gsea_score,ncol = 1)
    rownames(gsea_score) <- rownames(adj_mat)
    smoothed <- netSmooth(gsea_score, adj_mat, alpha = a)
    return(smoothed)
  }else if(inherits(seu, "VoltRon")){
    if (!("gsea_rat_norm" %in% colnames(Metadata(seu)))) {
      seu <- gseaCalc(seu, genes, assay)
    }
    adj <- as.matrix(as_adjacency_matrix(seu@graph$radius))
    gsea_score <- Metadata(seu)$gsea_rat_norm
    gsea_score <- as.matrix(gsea_score,ncol = 1)
    rownames(gsea_score) <- rownames(Metadata(seu))
    smoothed <- netSmooth(gsea_score, adj, alpha = a)
    return(smoothed)
  }else{
    stop("Input must be either of type Seurat or VoltRon")
  }
}  
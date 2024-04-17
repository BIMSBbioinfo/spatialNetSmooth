#' perform spatial smoothing
#' @param seu    Seurat object
#' @param genes csv-file indicating which genes to use (optional)
#' @param assay   which assay to use (optional)
#' @param a alpha for smoothing (optional)
#' @return vector with smoothed scores
#' @export
#' @importFrom netSmooth netSmooth
#' @import Seurat
#' @import dplyr

spatial_smooth <- function(seu, genes = "data/Datasets - Ikarus - Gene_lists.csv", assay = "Spatial", a = 0.8){
  seu <- gseaCalc(seu, genes, assay)
  adj_mat <- adj_matrix(seu)
  gsea_score <- seu@meta.data$gsea_rat_norm
  gsea_score <- matrix(gsea_score,ncol = 1)
  rownames(gsea_score) <- rownames(adj_mat)
  smoothed <- netSmooth(gsea_score, adj_mat, alpha = a)
  return(smoothed)
}

#' union smoothing
#'
#' @param seu Seurat-Object
#' @param genes csv-file indicating which genes to use (optional)
#' @param assay which assay to use (optional)
#' @param a alpha for smoothing (optional)
#' @param graph which neighbour-graph to use, "nn" or "snn"
#'
#' @return vector with smoothed scores
#' @export
#'
#' @examples union_smooth(seu, a = 0.4)
#' @importFrom netSmooth netSmooth
#' @import Seurat
#' @import dplyr
union_smooth <- function(seu, genes = "V1_Breast_Cancer_Block_A_Section_2_spatial/Datasets - Ikarus - Gene_lists.csv", assay = "Spatial", a = 0.8, graph = "nn"){
  seu <- gseaCalc(seu, genes, assay)
  if (graph== "snn")  neighbours <- as.matrix(seu@graphs$Spatial_snn)
  else neighbours <- as.matrix(seu@graphs$Spatial_nn)
  gsea_score <- seu@meta.data$gsea_rat_norm
  gsea_score <- matrix(gsea_score,ncol = 1)
  rownames(gsea_score) <- rownames(neighbours)
  adj_mat <- adj_matrix(seu)
  union <- (adj_mat | neighbours)*1
  gsea_smoothed <- netSmooth(gsea_score, union, alpha = a)
  return(gsea_smoothed)

}

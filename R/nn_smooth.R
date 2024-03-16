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


nn_smooth <- function(seu, genes = "V1_Breast_Cancer_Block_A_Section_2_spatial/Datasets - Ikarus - Gene_lists.csv", assay = "Spatial", a = 0.8, graph = "nn"){
  seu <- gseaCalc(seu, genes, assay)
  if (graph== "snn")  neighbours <- as.matrix(seu@graphs$Spatial_snn)
  else neighbours <- as.matrix(seu@graphs$Spatial_nn)
  gsea_score <- seu@meta.data$gsea_rat_norm
  gsea_score <- matrix(gsea_score,ncol = 1)
  rownames(gsea_score) <- rownames(neighbours)
  smoothed <- netSmooth(gsea_score, neighbours, alpha = a)
  return(smoothed)
}

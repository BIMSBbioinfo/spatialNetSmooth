#' nn_spatial_smooth
#'
#' @param seu Seurat-Object
#' @param genes csv-file indicating which genes to use (optional)
#' @param a1 alpha for nn/snn smoothing (optional)
#' @param a2 alpha for spatial smoothing (optional)
#' @param assay which assay to use (optional)
#' @param graph which neighbour-graph to use, "nn" or "snn"
#'
#' @return vector with smoothed score
#' @export
#'
#' @examples nn_spatial_smooth(seu, a1 = 0.5, a2 = 0.8, graph = "nn")
#' @importFrom netSmooth netSmooth
#' @import Seurat
#' @import dplyr

nn_spatial_smooth <- function(se, genes = "Datasets - Ikarus - Gene_lists.csv", assay = "Spatial", a1 = 0.8, a2 = 0.8, graph = "nn"){
  if(inherits(se, "Seurat")){
  se <- gseaCalc(se, genes, assay)
  graphs <- paste(assay, "_", graph, sep="")
  neighbours <- as.matrix(se@graphs[[graphs]])
  gsea_score <- se@meta.data$gsea_rat_norm
  gsea_score <- matrix(gsea_score,ncol = 1)
  rownames(gsea_score) <- rownames(neighbours)
  smoothed <- netSmooth(gsea_score, neighbours, alpha = a1)
  adj_mat <- adj_matrix(se)
  smoothed <- netSmooth(smoothed, adj_mat, alpha = a2)
  return(smoothed)
  }
  else if(inherits(se, "VoltRon")){
    
  }else{
    stop("Input must be either of type Seurat or VoltRon")
  }
}

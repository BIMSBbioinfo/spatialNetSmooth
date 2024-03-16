#' Create Adjacency Matrix
#' @param seu Seurat-Object
#' @param k distance to neighbours
#' @return The Adjacency Matrix
#' @keywords internal
#' @importFrom parallelDist parDist

adj_matrix <- function(seu, k = 7) {
  coordinates <- as.matrix(GetTissueCoordinates(seu, scale = "lowres"))
  dist_m <- parDist(coordinates, method = "euclidean")
  dist_m <- as.matrix(dist_m)
  adj_matrix <- matrix(unlist(sapply(dist_m, function(x) ifelse(x <= k, 1, 0))),ncol = ncol(dist_m), byrow = TRUE)
  rownames(adj_matrix) <- rownames(dist_m)
  colnames(adj_matrix) <- colnames(dist_m)
  return(adj_matrix)
}

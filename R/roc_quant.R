#' roc_quant
#'
#' @param truth vector of truth, 0 indicating normal cells, 1 cells with condition
#' @param scores vector of calculated scores
#'
#' @return return best threshold and accuracy
#' @export
#' @importFrom pROC roc coords
#' @examples roc_quant(c(0,1,1), c(-0.03, 0.6, 0.87))
roc_quant <- function(truth, scores){
  roc_c <-roc(response= truth, predictor=scores)
  coords <- coords(roc_c, "best", best.method="closest.topleft")
  best_threshold <- coords$threshold
  pred <- ifelse(scores >=best_threshold, 1, 0)
  conf_mat <- table(truth, pred)
  accuracy <- sum(diag(conf_mat))/sum(conf_mat)
  result <- NULL
  result$threshold <- best_threshold
  result$accuracy <- accuracy
  return(result)
}
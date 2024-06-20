#' F1_quant
#'
#' @param gsea_score calculated score for each cell
#' @param tumor vector for all cells with 1 for tumor and 0 for normal cells
#'
#' @return vector of F1-scores for different thresholds
#' @export
#' @importFrom caret confusionMatrix
#' @examples scores <-c(0.64, 0.048, 0.01)
#' tumor <- c(1, 0, 0)
#' F1_quant(scores, tumor)
F1_quant <- function(gsea_score, tumor){
  quant <- quantile(gsea_score, probs = seq(0,1,by =0.05))
  quant_min <- quant[2]
  quant_max <- quant[20]
  for(i in seq_along(gsea_score)) {
    if(gsea_score[i] < quant_min){
      gsea_score[i] <- quant_min
    }else if(gsea_score[i] > quant_max){
      gsea_score[i] <- quant_max
    }
  }
  vec <-NULL
  threshold <- quantile(gsea_score, probs = seq(0,1,by =0.1))
  for (g in threshold) {
    for(i in seq_along(gsea_score)) {
      if(gsea_score[i] < g){
        vec <- append(vec,0)
      }else if(gsea_score[i] > g){
        vec <- append(vec,1)
      }

    }
  }
  F1 <- NULL
  F1$threshold <- threshold
  n <- length(gsea_score)
  quant1 <- vec[1:n]
  quant2 <- vec[(n+1):(n*2)]
  quant3 <- vec[(n*2+1):(n*3)]
  quant4 <- vec[(n*3+1):(n*4)]
  quant5 <- vec[(n*4+1):(n*5)]
  quant6 <- vec[(n*5+1):(n*6)]
  quant7 <- vec[(n*6+1):(n*7)]
  quant8 <- vec[(n*7+1):(n*8)]
  quant9 <- vec[(n*8+1):(n*9)]
  quant10 <- vec[(n*9+1):(n*10)]
  quant11 <- vec[(n*10+1):(n*11)]

  tumor <- as.factor(tumor)
  quant1 <- as.factor(quant1)
  quant2 <- as.factor(quant2)
  quant3 <- as.factor(quant3)
  quant4 <- as.factor(quant4)
  quant5 <- as.factor(quant5)
  quant6 <- as.factor(quant6)
  quant7 <- as.factor(quant7)
  quant8 <- as.factor(quant8)
  quant9 <- as.factor(quant9)
  quant10 <- as.factor(quant10)
  quant11 <- as.factor(quant11)
  confusion_matrix <- confusionMatrix(quant1, tumor, positive="1")
  F1$score <- append(F1,confusion_matrix$byClass["F1"])
  confusion_matrix <- confusionMatrix(quant2, tumor, positive="1")
  F1$score <- append(F1,confusion_matrix$byClass["F1"])
  confusion_matrix <- confusionMatrix(quant3, tumor, positive="1")
  F1$score <- append(F1,confusion_matrix$byClass["F1"])
  confusion_matrix <- confusionMatrix(quant4, tumor, positive="1")
  F1$score <- append(F1,confusion_matrix$byClass["F1"])
  confusion_matrix <- confusionMatrix(quant5, tumor, positive="1")
  F1$score <- append(F1,confusion_matrix$byClass["F1"])
  confusion_matrix <- confusionMatrix(quant6, tumor, positive="1")
  F1$score <- append(F1,confusion_matrix$byClass["F1"])
  confusion_matrix <- confusionMatrix(quant7, tumor, positive="1")
  F1$score <- append(F1,confusion_matrix$byClass["F1"])
  confusion_matrix <- confusionMatrix(quant8, tumor, positive="1")
  F1$score <- append(F1,confusion_matrix$byClass["F1"])
  confusion_matrix <- confusionMatrix(quant9, tumor, positive="1")
  F1$score <- append(F1,confusion_matrix$byClass["F1"])
  confusion_matrix <- confusionMatrix(quant10, tumor, positive="1")
  F1$score <- append(F1,confusion_matrix$byClass["F1"])
  confusion_matrix <- confusionMatrix(quant11, tumor, positive="1")
  F1$score <- append(F1,confusion_matrix$byClass["F1"])
  return(F1)
}

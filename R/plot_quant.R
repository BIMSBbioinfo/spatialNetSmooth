#' plot_quant
#'
#' @param gsea_score
#' @param truth
#' @param threshold which percentile to take as threshold (1-10)
#' @param coordinates
#'
#' @return plot of truth of prediction
#' @export
#' @import ggplot2
#' @examples
plot_quant <- function(gsea_score, coordinates, truth, threshold){
  quant <- quantile(gsea_score, probs = seq(0,1,by =0.05))
  quant_min <- quant[2]
  quant_max <- quant[20]
  for(i in seq_along(gsea)) {
    if(gsea[i] < quant_min){
      gsea[i] <- quant_min
    }else if(gsea[i] > quant_max){
      gsea[i] <- quant_max
    }
  }

  thresh <- quantile(gsea_score, probs = seq(0,1,by =0.1))
  thresh <- thresh[threshold]
  for (i in seq_along(truth)) {
    if (gsea_score[i] >= thresh && tumor[i] == 1) {
      quant <- c(quant, 4)
    } else if(gsea_score[i] < thresh && tumor[i] == 0){
      quant <- c(quant, 1)
    }else if(gsea_score[i] < thresh && tumor[i] == 1){
      quant <- c(quant, 3)
    }else{
      quant <- c(quant, 2)
    }
  }
  # Create a data frame from the coordinates and values
  data_df <- data.frame(x = coordinates[, 1], y = coordinates[, 2], value = quant)

  # Simplify color mapping
  colors <- c("blue1", "lightblue", "orangered", "orangered4")  # Order corresponds to 0 and 1

  data_df$value <- as.numeric(data_df$value)
  # Plot using ggplot2
  p1 <- ggplot(data_df, aes(y, -x, color = factor(value))) +
    geom_point() +
    scale_color_manual(values = colors, labels = c("1" = "normal", "2" = "false tumor", "3" = "false normal", "4" = "tumor"))+
    ggtitle("truth of prediction")

  p1

}

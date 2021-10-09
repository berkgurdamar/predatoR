#' Impact Prediction
#'
#' @param final_df data.frame contains all the mutation information and the required information for impact prediction
#'
#' @import adabag
#' @return matrix
#' @export
#'

imp_prediction <- function(final_df){

  pred <- stats::predict(model, final_df)

  res <- pred$class

  probs <- c()
  res_types <- c()
  for(i in 1:length(res)){
    if(res[i] == 0){
      probs <- c(probs, pred$prob[i,1])
      res_types <- c(res_types, "Silent")
    }
    else{
      probs <- c(probs, pred$prob[i,2])
      res_types <- c(res_types, "Disease Causing")
    }
  }

  pred_df <- data.frame(Prediction = res_types,
                        Probability = probs)

  return(cbind(final_df[,1:5], pred_df))
}

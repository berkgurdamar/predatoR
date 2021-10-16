#' Impact Prediction
#'
#' Function that makes the final prediction
#'
#' This function takes a data.frame which contains all the required fields for impact prediction and compute
#' the prediction based on pre-computed Adaboost model.
#'
#' @param final_df data.frame contains all the mutation information and the required information for impact prediction
#'
#@import adabag
#' @return matrix contains the prediction results
#' @export
#'

imp_prediction <- function(final_df){

  # colnames(final_df)[4:5] <- c("X.Orig.Amino.Acid", "X.Mutant.Amino.Acid")
  prob <- stats::predict(caret_adaboost, final_df, type = "prob")

  res <- stats::predict(caret_adaboost, final_df)

  probs <- c()
  res_types <- c()
  for(i in 1:length(res)){
    if(res[i] == 0){
      probs <- c(probs, prob[i,1])
      res_types <- c(res_types, "Silent")
    }
    else{
      probs <- c(probs, prob[i,2])
      res_types <- c(res_types, "Disease Causing")
    }
  }

  pred_df <- data.frame(Prediction = res_types,
                        Probability = probs)

  return(cbind(final_df[,1:5], pred_df))
}

#' Impact Prediction
#'
#' Function that makes the final prediction
#'
#' This function takes a data.frame which contains all the required fields for impact prediction and compute
#' the prediction based on pre-computed Adaboost model.
#'
#' @param final_df data.frame contains all the required information for impact prediction
#'
#' @return matrix contains the prediction results
#'
#' @export
#'

impact_prediction <- function(final_df){

  final_df$Orig_AA <- as.factor(final_df$Orig_AA)
  final_df$Mut_AA <- as.factor(final_df$Mut_AA)
  final_df$degree_z_score <- as.numeric(final_df$degree_z_score)
  final_df$eigen_z_score <- as.numeric(final_df$eigen_z_score)
  final_df$shortest_path_z <- as.numeric(final_df$shortest_path_z)
  final_df$betweenness_scores_z <- as.numeric(final_df$betweenness_scores_z)
  final_df$clique_z_score <- as.numeric(final_df$clique_z_score)
  final_df$pagerank_z_score <- as.numeric(final_df$pagerank_z_score)
  final_df$syn_z <- as.numeric(final_df$syn_z)
  final_df$mis_z <- as.numeric(final_df$mis_z)
  final_df$pLI <- as.numeric(final_df$pLI)
  final_df$blosum62_scores <- as.numeric(final_df$blosum62_scores)
  final_df$kegg_pathway_number <- as.numeric(final_df$kegg_pathway_number)
  final_df$genic_intolerance <- as.numeric(final_df$genic_intolerance)
  final_df$go_terms <- as.numeric(final_df$go_terms)
  final_df$disgenet <- as.numeric(final_df$disgenet)
  final_df$gene_essentiality <- as.numeric(final_df$gene_essentiality)
  final_df$gtex <- as.numeric(final_df$gtex)
  final_df$ref_asa <- as.numeric(final_df$ref_asa)
  final_df$mut_asa <- as.numeric(final_df$mut_asa)
  final_df$asa_diff <- as.numeric(final_df$asa_diff)
  final_df$ref_hyd <- as.numeric(final_df$ref_hyd)
  final_df$mut_hyd <- as.numeric(final_df$mut_hyd)
  final_df$hyd_diff <- as.numeric(final_df$hyd_diff)

  prob <- stats::predict(caret_adaboost, final_df, type = "prob")

  res <- as.factor(ifelse(stats::predict(caret_adaboost, final_df, type = "prob")[,2] > 0.5, "1", "0"))

  probs <- c()
  res_types <- c()
  for(i in 1:length(res)){
    if(res[i] == 0){
      probs <- c(probs, prob[i,1])
      res_types <- c(res_types, "Neutral")
    }
    else{
      probs <- c(probs, prob[i,2])
      res_types <- c(res_types, "Pathogenic")
    }
  }

  pred_df <- data.frame(Prediction = res_types,
                        Probability = probs)

  return(cbind(final_df[,1:6], pred_df))
}

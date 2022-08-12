#' Impact Prediction
#'
#' Function that makes the final prediction
#'
#' This function takes a data.frame which contains all the required features and compute
#' the impact of a mutation.
#'
#' @param final_df data.frame contains all the required features for impact prediction
#' @param distance_cutoff distance cutoff for setting edges (default = 7)
#' @param network_approach network building approach; "all" (default) for using all atoms or "ca" for using carbon alpha atoms only
#'
#' @return data.frame contains the prediction results
#'
#' @export
#'

impact_prediction <- function(final_df, distance_cutoff = 7, network_approach = "all"){

  if(!network_approach %in% c("all", "ca")){
    stop("Network approach needs to be 'all' or 'ca'")

  }
  if(distance_cutoff == 7){

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

    if(network_approach == "all"){

      prob <- caret::predict.train(caret_adaboost_all, final_df, type = "prob")

      res <- as.factor(ifelse(caret::predict.train(caret_adaboost_all, final_df, type = "prob")[,2] > 0.4880486, "1", "0"))

    }else{

      prob <- caret::predict.train(caret_adaboost_ca, final_df, type = "prob")

      res <- as.factor(ifelse(caret::predict.train(caret_adaboost_ca, final_df, type = "prob")[,2] > 0.5505579, "1", "0"))

    }
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

  }else{
    stop("Predictions can be made using 7 angstrom cutoffs")
  }
}

#' Eigen Centrality Score
#'
#' Calculate Eigen Centrality Z-Scores of input position
#'
#' This function calculates Eigen Centrality scores of all nodes in the network.
#' Calculates and returns the Z-Scores of the input positions.
#'
#' @param connections_df data.frame contains all the edges
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return Eigen Centrality Z-Score of input position
#' @export
#'


eigen_centrality_score <- function(connections_df, filtered_info_df){

  final_df <- c()
  for(i in 1:length(connections_df)){
    total_scores <- c()
    idx <- unique(connections_df[[i]][grep("_CA_", connections_df[[i]][,1]), 1])
    for(j in idx){
      eigen_idx <- connections_df[[i]][which(connections_df[[i]][,1] == j), 2]
      total_scores <- c(total_scores, nrow(connections_df[[i]][connections_df[[i]][,1] %in% eigen_idx,]))
    }

    mean_of_scores <- mean(total_scores)
    sd_of_scores <- stats::sd(total_scores)

    z_scores <- (total_scores - mean_of_scores) / sd_of_scores
    final_df <- rbind(final_df, cbind(idx, z_scores))
  }
  final_df <- as.data.frame(final_df)

  z_final_scores <- as.numeric(sapply(paste0(filtered_info_df$Position, "_CA_", filtered_info_df$Chain),
                                      function(x) final_df$z_scores[which(final_df$idx == x)]))

  message(crayon::white(paste0("Eigen Centrality Score:", "\t\t", "DONE")))

  return(z_final_scores)
}

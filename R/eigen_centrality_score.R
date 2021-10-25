#' Eigen Centrality Score
#'
#' Calculate Eigen Centrality Z-Scores of input position
#'
#' This function calculates Eigen Centrality scores of all nodes in the network.
#' Calculates and returns the Z-Scores of the input positions.
#'
#' @param edge_list list contains separate edge data.frames for each chain
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return Eigen Centrality Z-Score of input position
#'
#' @export
#'


eigen_centrality_score <- function(edge_list, filtered_info_df){

  colnames(filtered_info_df)[1:5] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA")

  if(length(unique(filtered_info_df$PDB_ID)) > 1){
    stop(paste0("filtered_info_df should contain only one PDB entries"))
  }

  final_df <- c()
  for(i in 1:length(edge_list)){
    total_scores <- c()
    idx <- unique(edge_list[[i]][grep("_CA_", edge_list[[i]][,1]), 1])
    for(j in idx){
      eigen_idx <- edge_list[[i]][which(edge_list[[i]][,1] == j), 2]
      total_scores <- c(total_scores, nrow(edge_list[[i]][edge_list[[i]][,1] %in% eigen_idx,]))
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

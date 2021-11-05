#' Eigen Centrality Score
#'
#' Calculate Eigen Centrality Z-Scores of input position
#'
#' This function calculates the total number of connections of nodes to which a node is connected.
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

    edge_list_filtered <- as.data.frame(edge_list[[i]])

    df.g <- igraph::graph.data.frame(d = edge_list_filtered, directed = FALSE)

    eigen_centrality_scores <- igraph::eigen_centrality(df.g)$vector

    eigen_z_score <- (eigen_centrality_scores - mean(eigen_centrality_scores)) / stats::sd(eigen_centrality_scores)

    final_df <- rbind(final_df, cbind(names(eigen_z_score), eigen_z_score))

  }
  final_df <- as.data.frame(final_df)
  colnames(final_df)[1] <- "res_name"

  z_final_scores <- as.numeric(sapply(paste0(filtered_info_df$Position, "_CA_", filtered_info_df$Chain),
                                      function(x) final_df$eigen_z_score[which(final_df$res_name == x)]))

  message(crayon::white(paste0("Eigen Centrality Score:", "\t\t", "DONE")))

  return(z_final_scores)
}

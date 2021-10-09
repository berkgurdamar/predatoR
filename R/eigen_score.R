#' Eigen Centrality Score
#'
#' @param connections_df data.frame contains all the edges
#' @param Position Mutation position
#' @param atom_matrix matrix that contains all the atoms in the PDB file
#'
#' @return numeric
#' @export
#'

eigen_score <- function(connections_df, Position, atom_matrix){
  idx <- strsplit(Position, "_")[[1]]

  eigen_idx <- connections_df$connections[which(connections_df$node_name == paste0(idx[1], "_CA_", idx[2]))]
  scores <- nrow(connections_df[connections_df$node_name %in% eigen_idx,])

  ca_df <- atom_matrix[atom_matrix$elety == "CA",]
  total_scores <- c()
  for(k in 1:nrow(ca_df)){
    eigen_idx <- connections_df$connections[which(connections_df$node_name == paste0(ca_df$resno[k], "_CA_", ca_df$chain[k]))]
    total_scores <- c(total_scores, nrow(connections_df[connections_df$node_name %in% eigen_idx,]))
  }
  mean_of_scores <- mean(total_scores)
  sd_of_scores <- stats::sd(total_scores)

  eigen_z_score <- (scores - mean_of_scores) / sd_of_scores

  return(eigen_z_score)
}

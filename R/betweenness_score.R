#' Betweenness Score
#'
#' Calculate Betweenness Z-Scores of input position
#'
#' This function calculates Betweenness scores of all nodes in the network via \code{betweenness} function
#' of igraph package. Calculates and returns the Z-Scores of the input positions.
#'
#' @param connections_df data.frame contains all the edges calculated by \code{PDB2connections} function
#' @param Position Mutation position
#'
#' @return Betweenness Z-Scores of the input position
#' @export
#'

betweenness_score <- function(connections_df, Position){

  idx <- strsplit(Position, "_")[[1]]

  connections_df <- as.data.frame(connections_df[[which(names(connections_df) == idx[2])]])

  df.g <- igraph::graph.data.frame(d = connections_df, directed = FALSE)

  all_betwenness <- igraph::betweenness(df.g, directed = F)
  mean_betwenness <- mean(all_betwenness)
  sd_betwenness <- stats::sd(all_betwenness)

  idx_betwenness <- all_betwenness[paste0(idx[1], "_CA_", idx[2])]

  betwenness_scores_z <- (idx_betwenness - mean_betwenness) / sd_betwenness

  return(betwenness_scores_z)
}

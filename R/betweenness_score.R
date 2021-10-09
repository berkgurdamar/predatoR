#' Betweenness Score
#'
#' @param connections_df data.frame contains all the edges
#' @param Position Mutation position
#'
#' @return numeric
#' @export
#'

betweenness_score <- function(connections_df, Position){

  idx <- strsplit(Position, "_")[[1]]

  df.g <- igraph::graph.data.frame(d = connections_df, directed = FALSE)

  all_betwenness <- igraph::betweenness(df.g, directed = F)
  mean_betwenness <- mean(all_betwenness)
  sd_betwenness <- stats::sd(all_betwenness)

  idx_betwenness <- all_betwenness[paste0(idx[1], "_CA_", idx[2])]

  betwenness_scores_z <- (idx_betwenness - mean_betwenness) / sd_betwenness

  return(betwenness_scores_z)
}

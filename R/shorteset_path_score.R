#' Shortest Path Score Function
#'
#' Calculates shortest path z-score of mutation index amino acid from list of edges
#'
#' @param connections_df data frame contains list of edges
#' @param Position Mutation position
#'
#' @return integer
#' @export
#'

shorteset_path_score <- function(connections_df, Position){

  idx <- strsplit(Position, "_")[[1]]

  df.g <- igraph::graph.data.frame(d = connections_df, directed = FALSE)

  distMatrix <- igraph::shortest.paths(df.g, v=igraph::V(df.g), to=igraph::V(df.g))
  distMatrix[which(distMatrix == "Inf")] <- 0

  mean_shortest_path_length <- mean(rowSums(distMatrix))
  sd_shortest_path_length <- stats::sd(rowSums(distMatrix))

  shortest_paths <- rowSums(distMatrix)[paste0(idx[1], "_CA_", idx[2])]

  shortest_path_z <- (shortest_paths - mean_shortest_path_length) / sd_shortest_path_length

  return(shortest_path_z)
}

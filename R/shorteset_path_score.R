#' Shortest Path Score Function
#'
#' Calculate Shortest Path Z-Scores of input position
#'
#' This function calculates Shortest Path lengths of all nodes in the network via \code{shortest.paths} function
#' of igraph package. Sum shortest path lengts of every nodes, calculates and returns the Z-Scores of the input positions.
#'
#' @param connections_df data frame contains list of edges
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return Shortest Path Z-Scores of input position
#' @export
#'


shorteset_path_score <- function(connections_df, filtered_info_df){

  final_df <- c()
  for(i in 1:length(connections_df)){
    connections_df_filtered <- as.data.frame(connections_df[[i]])
    df.g <- igraph::graph.data.frame(d = connections_df_filtered, directed = FALSE)

    distMatrix <- igraph::shortest.paths(df.g, v=igraph::V(df.g), to=igraph::V(df.g))
    distMatrix[which(distMatrix == "Inf")] <- 0

    mean_shortest_path_length <- mean(rowSums(distMatrix))
    sd_shortest_path_length <- stats::sd(rowSums(distMatrix))

    shortest_path_z <- (rowSums(distMatrix) - mean_shortest_path_length) / sd_shortest_path_length

    final_df <- rbind(final_df, cbind(names(shortest_path_z), shortest_path_z))
    }

  final_df <- as.data.frame(final_df)

  shortest_z_final <- as.numeric(sapply(paste0(filtered_info_df$Position, "_CA_", filtered_info_df$Chain),
                                        function(x) final_df$shortest_path_z[which(final_df$V1 == x)]))

  message(crayon::white(paste0("Number of Shortest Paths:", "\t", "DONE")))

  return(shortest_z_final)
}

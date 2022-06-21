#' Shortest Path Score
#'
#' Calculate Shortest Path Z-Scores of input position
#'
#' This function calculates Shortest Path lengths of all nodes in the network via \code{shortest.paths} function
#' of igraph package. Sum shortest path lengths of every nodes, calculates and returns the Z-Scores of the input positions.
#'
#' @param edge_list list contains separate edge data.frames for each chain
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return Shortest Path Z-Scores of input position
#'
#' @export
#'


shorteset_path_score <- function(edge_list, filtered_info_df){

  colnames(filtered_info_df)[1:5] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA")

  if(length(unique(filtered_info_df$PDB_ID)) > 1){
    stop(paste0("filtered_info_df should contain only one PDB entries"))
  }

  final_df <- c()
  for(i in 1:length(edge_list)){
    edge_list_filtered <- as.data.frame(edge_list[[i]])
    df.g <- igraph::graph.data.frame(d = edge_list_filtered, directed = FALSE)

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

  message(crayon::white(paste0("Shortest Path Score:", "\t\t", "DONE")))

  return(shortest_z_final)
}

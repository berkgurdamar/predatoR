#' Betweenness Score
#'
#' Calculate Betweenness Z-Scores of input position
#'
#' This function calculates Betweenness scores of all nodes in the network via \code{betweenness} function
#' of igraph package. Calculates and returns the Z-Scores of the input positions.
#'
#' @param connections_df data.frame contains all the edges calculated by \code{PDB2connections} function
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return Betweenness Z-Scores of the input position
#' @export
#'

betweenness_score <- function(connections_df, filtered_info_df){

  final_df <- c()
  for(i in 1:length(connections_df)){

    connections_df_filtered <- as.data.frame(connections_df[[i]])

    df.g <- igraph::graph.data.frame(d = connections_df_filtered, directed = FALSE)

    all_betwenness <- igraph::betweenness(df.g, directed = F)
    mean_betwenness <- mean(all_betwenness)
    sd_betwenness <- stats::sd(all_betwenness)

    final_betweenness <- (all_betwenness - mean_betwenness) / sd_betwenness

    final_df <- rbind(final_df, cbind(names(final_betweenness), final_betweenness))

  }

  final_df <- as.data.frame(final_df)

  betweenness_z_final <- as.numeric(sapply(paste0(filtered_info_df$Position, "_CA_", filtered_info_df$Chain),
                                           function(x) final_df$final_betweenness[which(final_df$V1 == x)]))

  message(crayon::white(paste0("Betweenness Score:", "\t\t", "DONE")))

  return(betweenness_z_final)
}


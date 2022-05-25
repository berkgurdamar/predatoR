#' Betweenness Score
#'
#' Calculate Betweenness Z-Scores of input position
#'
#' This function calculates Betweenness scores of all nodes in the network via \code{betweenness} function
#' of igraph package. Calculates and returns the Z-Scores of the input positions.
#'
#' @param edge_list list contains separate edge data.frames for each chain
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return Betweenness Z-Scores of the input position
#'
#' @export
#'

betweenness_score <- function(edge_list, filtered_info_df){

  colnames(filtered_info_df)[1:5] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA")

  if(length(unique(filtered_info_df$PDB_ID)) > 1){
    stop(paste0("filtered_info_df should contain only one PDB entries"))
  }


  final_df <- c()
  for(i in 1:length(edge_list)){

    edge_list_filtered <- as.data.frame(edge_list[[i]])

    df.g <- igraph::graph.data.frame(d = edge_list_filtered, directed = FALSE)

    all_betwenness <- igraph::betweenness(df.g, directed = F)
    all_betwenness[which(all_betwenness == "Inf")] <- 0
    all_betwenness[which(all_betwenness == "NaN")] <- 0

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

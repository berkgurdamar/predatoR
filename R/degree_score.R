#' Degree Score
#'
#' Calculates and returns the Degree Z-Scores of the input positions.
#'
#' @param edge_list list contains separate edge data.frames for each chain
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return Degree Z-Scores of input position
#' @export
#'


degree_score <- function(edge_list, filtered_info_df){

  colnames(filtered_info_df)[1:5] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA")

  if(length(unique(filtered_info_df$PDB_ID)) > 1){
    stop(paste0("filtered_info_df should contain only one PDB entries"))
  }

  final_df <- c()
  for(i in 1:length(edge_list)){

    edge_list_filtered <- as.data.frame(edge_list[[i]])

    df.g <- igraph::graph.data.frame(d = edge_list_filtered, directed = FALSE)

    degree_scores <- igraph::degree(df.g,
                                    v = igraph::V(df.g))

    degree_z_scores <- (degree_scores - mean(degree_scores)) / stats::sd(degree_scores)

    final_df <- rbind(final_df, cbind(names(degree_z_scores), degree_z_scores))

  }

  final_df <- as.data.frame(final_df)
  colnames(final_df)[1] <- "res_name"

  z_final_scores <- as.numeric(sapply(paste0(filtered_info_df$Position, "_CA_", filtered_info_df$Chain),
                                      function(x) final_df$degree_z_scores[which(final_df$res_name == x)]))

  message(crayon::white(paste0("Degree Score:", "\t\t\t", "DONE")))

  return(z_final_scores)

}

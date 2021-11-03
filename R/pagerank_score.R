#' PageRank Score
#'
#' \code{pagerank_score} function assigns an importance score to every node via
#' \code{page_rank} function of igraph package..
#'
#' @param edge_list list contains separate edge data.frames for each chain
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return PageRank Z-Scores of input position
#' @export
#'


pagerank_score <- function(edge_list, filtered_info_df){

  colnames(filtered_info_df)[1:5] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA")

  if(length(unique(filtered_info_df$PDB_ID)) > 1){
    stop(paste0("filtered_info_df should contain only one PDB entries"))
  }

  final_df <- c()
  for(i in 1:length(edge_list)){

    edge_list_filtered <- as.data.frame(edge_list[[i]])

    df.g <- igraph::graph.data.frame(d = edge_list_filtered, directed = FALSE)

    pagerank_scores <- igraph::page_rank(df.g)$vector

    pagerank_z_scores <- (pagerank_scores - mean(pagerank_scores)) / stats::sd(pagerank_scores)

    final_df <- rbind(final_df, cbind(names(pagerank_z_scores), pagerank_z_scores))

  }

  final_df <- as.data.frame(final_df)
  colnames(final_df)[1] <- "res_name"

  z_final_scores <- as.numeric(sapply(paste0(filtered_info_df$Position, "_CA_", filtered_info_df$Chain),
                                      function(x) final_df$pagerank_z_scores[which(final_df$res_name == x)]))

  message(crayon::white(paste0("PageRank Score:", "\t\t\t", "DONE")))

  return(z_final_scores)
}

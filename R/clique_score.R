#' Clique Score
#'
#' Calculates and returns the Clique Z-Scores of the input positions.
#'
#' @param edge_list list contains separate edge data.frames for each chain
#' @param filtered_info_df input data.frame which contain only one PDB entries
#' @param n_threads number of threads (default = NULL)
#' @param single_run should be set as TRUE when using \code{clique_score} function alone (default = TRUE)
#'
#' @return Clique Z-Scores of input position
#' @export
#'

clique_score <- function(edge_list, filtered_info_df, n_threads = NULL, single_run = TRUE){

  colnames(filtered_info_df)[1:5] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA")

  if(single_run == TRUE){

    close_parallel <- TRUE

    if(is.null(n_threads)){

      n.cores <- parallel::detectCores() - 1

    }else{

      n.cores <- n_threads

    }

    if(length(unique(filtered_info_df$PDB_ID)) > 1){
      stop(paste0("filtered_info_df should contain only one PDB entries"))
    }

    clusters <- parallel::makeCluster(n.cores, strategy = "sequential")
    doParallel::registerDoParallel(clusters)

  }else{

    close_parallel <- FALSE

    if(length(unique(filtered_info_df$PDB_ID)) > 1){
      stop(paste0("filtered_info_df should contain only one PDB entries"))
    }
  }

  `%dopar%` <- foreach::`%dopar%`

  j <- ""
  final_df <- c()
  for(i in 1:length(edge_list)){

    edge_list_filtered <- as.data.frame(edge_list[[i]])

    all_perm <- choose(length(unique(edge_list_filtered$node_name)), 2) * factorial(2)

    clique_score <- foreach::foreach(j = 1:length(unique(edge_list_filtered$node_name)), .combine = c) %dopar% {

      idx <- edge_list_filtered[which(edge_list_filtered[,1] == unique(edge_list_filtered$node_name)[j]), 2]
      cliques <- edge_list_filtered[intersect(which(edge_list_filtered$node_name %in% idx), which(edge_list_filtered$connections %in% idx)),]

      nrow(cliques)/all_perm
    }

    clique_z_score <- (clique_score - mean(clique_score)) / stats::sd(clique_score)
    final_df <- rbind(final_df, cbind(unique(edge_list_filtered$node_name), clique_z_score))

  }

  final_df <- as.data.frame(final_df)
  colnames(final_df)[1] <- "res_name"

  z_final_scores <- as.numeric(sapply(paste0(filtered_info_df$Position, "_CA_", filtered_info_df$Chain),
                                      function(x) final_df$clique_z_score[which(final_df$res_name == x)]))

  if(close_parallel == TRUE){

    parallel::stopCluster(clusters)

  }

  message(crayon::white(paste0("Clique Score:", "\t\t\t", "DONE")))

  return(z_final_scores)

}


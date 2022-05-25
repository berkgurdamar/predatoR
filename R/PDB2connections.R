#' Creating edge list from atom matrix
#'
#' Creates all the edges from the input atom matrix created by \code{PDB_read} .
#'
#' This function calculates the distances between all the atoms in the atom matrix and create edges
#'  if the distance between two atoms is less or equal to 7 Angstrom.
#'
#' @param atom_matrix matrix that contains all the atoms in the PDB file
#' @param filtered_info_df input data.frame which contain only one PDB entries
#' @param n_threads number of threads (default = NULL)
#' @param single_run should be set as TRUE when using \code{PDB2connections} function alone (default = TRUE)
#'
#' @return list contains separate edge data.frames for each chain
#'
#' @export
#'

PDB2connections <- function(atom_matrix, filtered_info_df, n_threads = NULL, single_run = TRUE){

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

<<<<<<< HEAD
=======

>>>>>>> ee79646b65c39e937123c1002f0ae00cbbfee369
  `%dopar%` <- foreach::`%dopar%`

  j <- ""
  edge_list <- list()
  for(index in 1:length(unique(filtered_info_df$Chain))){
    filtered_atom_matrix <- atom_matrix[atom_matrix$chain == unique(filtered_info_df$Chain)[index],]

    chain_df <- foreach::foreach(j = 1:nrow(filtered_atom_matrix), .combine = rbind) %dopar% {

      node_name <- c()
      connections <- c()

      others <- filtered_atom_matrix[c("x", "y", "z")][-j,]
      idx <- matrix(as.numeric(filtered_atom_matrix[c("x", "y", "z")][j,]), nrow(others), 3, byrow = T)

      dist_mat <- sqrt(rowSums((others - idx)^2))
      names(dist_mat) <- paste0(filtered_atom_matrix[names(dist_mat), "resno"],
                                "_", filtered_atom_matrix[names(dist_mat), "elety"],
                                "_", filtered_atom_matrix[names(dist_mat), "chain"])

      node_name <- c(node_name, as.character(rep(paste0(filtered_atom_matrix$resno[j],
                                                        "_", filtered_atom_matrix$elety[j],
                                                        "_", filtered_atom_matrix$chain[j]),
                                                 length(names(dist_mat[which(dist_mat <= 7)])))))

      connections <- c(connections, names(dist_mat[which(dist_mat <= 7)]))

      cbind(node_name, connections)

    }
    edge_list[[length(edge_list) + 1]] <- chain_df

  }

  names(edge_list) <- unique(filtered_info_df$Chain)

  if(close_parallel == TRUE){

    parallel::stopCluster(clusters)

  }
  message(crayon::white(paste0("Edge List:", "\t\t\t", "DONE")))

  return(edge_list)
}

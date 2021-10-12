#' Creating edge list from atom matrix
#'
#' Creates all the edges from the input atom matrix created by \code{PDB_read} function
#'
#' This function calculates the distances between all the atoms in the atom matrix and create edges
#'  if the distance between two atoms is less or equal to 7 Angstrom.
#'
#' @param atom_matrix matrix that contains all the atoms in the PDB file
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return list containing all the edges of each chain
#' @export
#'

PDB2connections <- function(atom_matrix, filtered_info_df){

  # doFuture::registerDoFuture()
  # n.cores <- parallel::detectCores() - 1
  #
  # my.cluster <- parallel::makeCluster(n.cores)
  # future::plan(future::cluster, workers = my.cluster)
  #
  `%dopar%` <- foreach::`%dopar%`

  j <- ""
  connections_df <- list()
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
    connections_df[[length(connections_df) + 1]] <- chain_df

  }

  names(connections_df) <- unique(filtered_info_df$Chain)

  # parallel::stopCluster(my.cluster)

  message(crayon::white(paste0("List of Edges:", "\t\t\t", "DONE")))

  return(connections_df)
}

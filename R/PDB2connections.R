#' Creating edge list from atom matrix
#'
#' Creates all the edges from the input atom matrix created by \code{PDB_read} function
#'
#' This function calculates the distances between all the atoms in the atom matrix and create edges
#'  if the distance between two atoms is less or equal to 7 Angstrom.
#'
#' @param atom_matrix matrix that contains all the atoms in the PDB file
#'
#' @return list containing all the edges of each chain
#' @export
#'

PDB2connections <- function(atom_matrix){


  doFuture::registerDoFuture()
  n.cores <- parallel::detectCores() - 1

  my.cluster <- parallel::makeCluster(n.cores)
  future::plan(future::cluster, workers = my.cluster)

  `%dopar%` <- foreach::`%dopar%`

  index = ""

  connections_df <- foreach::foreach(index = 1:length(unique(atom_matrix$chain))) %dopar% {
    filtered_atom_matrix <- atom_matrix[atom_matrix$chain == unique(atom_matrix$chain)[index],]
    node_name <- c()
    connections <- c()
    for(i in 1:nrow(filtered_atom_matrix)){
      others <- filtered_atom_matrix[c("x", "y", "z")][-i,]
      idx <- matrix(as.numeric(filtered_atom_matrix[c("x", "y", "z")][i,]), nrow(others), 3, byrow = T)

      dist_mat <- sqrt(rowSums((others - idx)^2))
      names(dist_mat) <- paste0(filtered_atom_matrix[names(dist_mat), "resno"],
                                "_", filtered_atom_matrix[names(dist_mat), "elety"],
                                "_", filtered_atom_matrix[names(dist_mat), "chain"])

      node_name <- c(node_name, as.character(rep(paste0(filtered_atom_matrix$resno[i],
                                                        "_", filtered_atom_matrix$elety[i],
                                                        "_", filtered_atom_matrix$chain[i]),
                                                 length(names(dist_mat[which(dist_mat <= 7)])))))

      connections <- c(connections, names(dist_mat[which(dist_mat <= 7)]))

    }
    cbind(node_name, connections)
  }

  names(connections_df) <- unique(atom_matrix$chain)
  # connections_df <- foreach::foreach(index = 1:nrow(atom_matrix), .combine = rbind) %dopar% {
  #
  #   node_name <- c()
  #   connections <- c()
  #
  #   others <- atom_matrix[c("x", "y", "z")][-index,]
  #   idx <- matrix(as.numeric(atom_matrix[c("x", "y", "z")][index,]), nrow(others), 3, byrow = T)
  #
  #   dist_mat <- sqrt(rowSums((others - idx)^2))
  #   names(dist_mat) <- paste0(atom_matrix[names(dist_mat), "resno"],
  #                             "_", atom_matrix[names(dist_mat), "elety"],
  #                             "_", atom_matrix[names(dist_mat), "chain"])
  #
  #   node_name <- c(node_name, as.character(rep(paste0(atom_matrix$resno[index],
  #                                                     "_", atom_matrix$elety[index],
  #                                                     "_", atom_matrix$chain[index]),
  #                                              length(names(dist_mat[which(dist_mat <= 7)])))))
  #
  #   connections <- c(connections, names(dist_mat[which(dist_mat <= 7)]))
  #
  #   cbind(node_name, connections)
  #
  # }

  parallel::stopCluster(my.cluster)

  message(crayon::white(paste0("List of Edges:", "\t\t\t", "DONE")))

  return(connections_df)
}

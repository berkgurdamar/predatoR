#' Creating edge list from atom matrix
#'
#' @param atom_matrix matrix that contains all the atoms in the PDB file
#'
#' @return data.frame
#' @export
#'

PDB2connections <- function(atom_matrix){


  doFuture::registerDoFuture()
  n.cores <- parallel::detectCores() - 1

  my.cluster <- parallel::makeCluster(n.cores)
  future::plan(future::cluster, workers = my.cluster)

  `%dopar%` <- foreach::`%dopar%`

  index = 1:nrow(atom_matrix)
  connections_df <- foreach::foreach(index = 1:nrow(atom_matrix), .combine = rbind) %dopar% {

    node_name <- c()
    connections <- c()

    others <- atom_matrix[c("x", "y", "z")][-index,]
    idx <- matrix(as.numeric(atom_matrix[c("x", "y", "z")][index,]), nrow(others), 3, byrow = T)

    dist_mat <- sqrt(rowSums((others - idx)^2))
    names(dist_mat) <- paste0(atom_matrix[names(dist_mat), "resno"],
                              "_", atom_matrix[names(dist_mat), "elety"],
                              "_", atom_matrix[names(dist_mat), "chain"])

    node_name <- c(node_name, as.character(rep(paste0(atom_matrix$resno[index],
                                                      "_", atom_matrix$elety[index],
                                                      "_", atom_matrix$chain[index]),
                                               length(names(dist_mat[which(dist_mat <= 7)])))))

    connections <- c(connections, names(dist_mat[which(dist_mat <= 7)]))

    cbind(node_name, connections)

  }

  parallel::stopCluster(my.cluster)

  message(crayon::white(paste0("Creating Atom Matrix:", "\t\t", "DONE")))

  return(as.data.frame(connections_df))
}

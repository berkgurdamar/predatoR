#' Mutation impact Prediction
#'
#' Main function for mutation impact prediction
#'
#' With an input data.frame which contains 'PDB_ID', 'Chain', 'Position',
#' 'Reference Amino Acid', 'Mutated Amino Acid' and 'Gene_Name' (optional) information respectively,
#' make a prediction about the mutation by using pre computed adaboost model
#' and classifies the mutation as 'Disease Causing' or 'Silent'.
#'
#' @param info_df data.frame containing 'PDB_ID', 'Chain', 'Position',
#' 'Reference Amino Acid', 'Mutated Amino Acid' and 'Gene_Name' (optional) information respectively.
#' @param PDB_path PDB file path (default = NULL)
#' @param n_threads number of threads (default = NULL)
#' @param gene_name_info whether there is gene name information in the input or not (default = TRUE)
#'
#' @return data.frame which contains prediction results
#' @export
#'

predatoR <- function(info_df, PDB_path = NULL, n_threads = NULL, gene_name_info = TRUE){

  if(is.data.frame(info_df) == F){
    stop("Input should be a data.frame.")
  }
  else if(gene_name_info == TRUE){
    if(ncol(info_df) != 6){
      stop("Input data.frame should contain 6 columns; PDB_ID, Chain, Position, Orig_AA and Mut_AA, Gene_Name respectively.")
    }
    else{
      colnames(info_df) <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA", "Gene_Name")
      info_df$Orig_AA <- toupper(info_df$Orig_AA)
      info_df$Mut_AA <- toupper(info_df$Mut_AA)
    }
  }
  else if(gene_name_info == FALSE){
    if(ncol(info_df) != 5){
      stop("Input data.frame should contain 5 columns; PDB_ID, Chain, Position, Orig_AA and Mut_AA respectively.")
    }
    else{
      colnames(info_df) <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA")
      info_df$Orig_AA <- toupper(info_df$Orig_AA)
      info_df$Mut_AA <- toupper(info_df$Mut_AA)
    }
  }
  else if(length(setdiff(unique(unlist(c(info_df[4], info_df[5]))), colnames(blosum_data))) > 0){

    false_name <- setdiff(unique(unlist(c(info_df[4], info_df[5]))), colnames(blosum_data))

    stop(paste0(paste0(false_name, collapse = ", "),
                " couldn't find in the amino acid names (Amino acid names should be 3 letter codes)"))
  }

  doFuture::registerDoFuture()

  if(is.null(n_threads) == TRUE){
    n.cores <- parallel::detectCores() - 1
  }else{
    n.cores <- n_threads
  }

  my.cluster <- parallel::makeCluster(n.cores)
  future::plan(future::cluster, workers = my.cluster)

  final_df <- c()
  for(i in unique(info_df$PDB_ID)){

    filtered_info_df <- info_df[info_df$PDB_ID == i,]

    atom_matrix <- read_PDB(i, PDB_path = PDB_path)

    if(is.data.frame(atom_matrix) == FALSE){
      next
    }
    message(crayon::white(paste0("PDB ID:", "\t\t\t\t", i, "\n",
                                 "Position(s):", "\t\t\t", paste0(filtered_info_df$Position, collapse = ", "))))

    removed_idx <- c()
    for(j in 1:nrow(filtered_info_df)){
      if(!(filtered_info_df$Position[j] %in% unique(atom_matrix$resno))){
        message(crayon::white(paste0("Residue ", filtered_info_df$Position[j],
                                     " is not included in the PDB structure, it will be removed from the query")))
        removed_idx <- c(removed_idx, j)
      }
      else if(length(unique(atom_matrix$resid[atom_matrix$resno == filtered_info_df$Position[j] & atom_matrix$chain == filtered_info_df$Chain[j]])) != 1){
        message(crayon::white(paste0("There are multiple amino acids for residue ", filtered_info_df$Position[j],
                                     " in the PDB file, it will be removed from the query")))
        removed_idx <- c(removed_idx, j)
        }
      else if(unique(atom_matrix$resid[atom_matrix$resno == filtered_info_df$Position[j] & atom_matrix$chain == filtered_info_df$Chain[j]]) != filtered_info_df$Orig_AA[j]){
        message(crayon::white(paste0("Residue ", filtered_info_df$Position[j], " is not ",
                                     filtered_info_df$Orig_AA[j], " in the PDB structure, it will be removed from the query")))
        removed_idx <- c(removed_idx, j)
      }
    }

    if(length(removed_idx) > 0){
    filtered_info_df <- filtered_info_df[-removed_idx,]
      if(nrow(filtered_info_df) == 0){
        next
      }
    }

    connections_df <- PDB2connections(atom_matrix, filtered_info_df)

    ### eigen centrality

    filtered_info_df$eigen_z_score <- eigen_centrality_score(connections_df, filtered_info_df)

    ### total shortest paths

    filtered_info_df$shortest_path_z <- shorteset_path_score(connections_df, filtered_info_df)

    ### betweenness scores

    filtered_info_df$betwenness_scores_z <- betweenness_score(connections_df, filtered_info_df)

    ### gnomad

    gnomad_result <- gnomad_scores(i, filtered_info_df)
    filtered_info_df <- gnomad_result[[1]]
    gene_name <- gnomad_result[[2]]

    if(gene_name == "no_name"){
      next
    }

    ### BLOSUM62

    filtered_info_df$blosum62_scores <- as.numeric(BLOSUM62_score(filtered_info_df))

    ### kegg pathway

    filtered_info_df$kegg_pathway_number <- as.numeric(KEGG_pathway_number(gene_name))

    ### genic intolerance

    filtered_info_df$genic_intolerance <- as.numeric(genic_intolerance(gene_name))

    final_df <- rbind(final_df, filtered_info_df)
  }

  ### prediction

  final_df <- stats::na.omit(final_df)

  if(nrow(final_df) == 0){

    stop("There is no input for prediction after NA omit")

  }else{

  prediction_result <- impact_prediction(final_df)

  }

  parallel::stopCluster(my.cluster)

  return(prediction_result)

}

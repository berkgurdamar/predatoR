#' Mutation impact Prediction
#'
#' Main function for mutation impact prediction
#'
#' With an input data.frame which contains 'PDB_ID', 'Chain', 'Position',
#' 'Reference Amino Acid' and 'Mutated Amino Acid' information respectively,
#' make a prediction about the mutation by using pre computed adaboost model
#' and classifies the mutation as 'Disease Causing' or 'Silent'.
#'
#' @param info_df data.frame containing 'PDB_ID', 'Chain', 'Position',
#' 'Reference Amino Acid' and 'Mutated Amino Acid' information respectively.
#'
#' @return data.frame which contains prediction results
#' @export
#'

PredImption <- function(info_df){

  colnames(info_df) <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA")

  final_df <- c()
  for(i in unique(info_df$PDB_ID)){

    filtered_info_df <- info_df[info_df$PDB_ID == i,]

    atom_matrix <- PDB_read(i)

    message(crayon::white(paste0("PDB ID:", "\t\t\t\t", i, "\n",
                                 "Position(s):", "\t\t\t", paste0(filtered_info_df$Position, collapse = ", "))))

    for(j in 1:nrow(filtered_info_df)){
      if(unique(atom_matrix$resid[atom_matrix$resno == filtered_info_df$Position[j] & atom_matrix$chain == filtered_info_df$Chain[j]]) != filtered_info_df$Orig_AA[j]){
        stop(paste0("Residue ", filtered_info_df$Position[j], " is not ", filtered_info_df$Orig_AA[j], " in the PDB structure."))
      }
      if(!(filtered_info_df$Position[j] %in% unique(atom_matrix$resno))){
        stop(paste0("Residue ", filtered_info_df$Position[j], " is not included in the PDB structure."))
      }
    }

    connections_df <- PDB2connections(atom_matrix)

    ### eigen centrality

    filtered_info_df$eigen_z_score <- sapply(paste0(filtered_info_df$Position, "_",filtered_info_df$Chain), function(x) eigen_score(connections_df, x, atom_matrix))

    message(crayon::white(paste0("Eigen Centrality Score:", "\t\t", "DONE")))

    ### total shortest paths

    filtered_info_df$shortest_path_z <- sapply(paste0(filtered_info_df$Position, "_",filtered_info_df$Chain), function(x) shorteset_path_score(connections_df, x))

    message(crayon::white(paste0("Number of Shortest Paths:", "\t", "DONE")))

    ### betweenness scores

    filtered_info_df$betwenness_scores_z <- sapply(paste0(filtered_info_df$Position, "_",filtered_info_df$Chain), function(x) betweenness_score(connections_df, x))

    message(crayon::white(paste0("Betweenness Score:", "\t\t", "DONE")))

    ### mutation type

    filtered_info_df$mut_type <- as.factor(mutation_type(filtered_info_df))

    ### gnomad

    gnomad_res <- gnomad(i, filtered_info_df)
    filtered_info_df <- gnomad_res[[1]]
    gene_name <- gnomad_res[[2]]

    ### BLOSUM62

    filtered_info_df$blosum62_scores <- blosum62_score(filtered_info_df)

    ### kegg pathway

    filtered_info_df$kegg_pathway_number <- as.numeric(kegg_pathway_number(gene_name))

    ### genic intolerance

    filtered_info_df$genic_intolerance <- as.numeric(genic_intolerance(gene_name))

    # close(pb)
    final_df <- rbind(final_df, filtered_info_df)
  }
  ### prediction
  final_df <- stats::na.omit(final_df)

  if(nrow(final_df) == 0){

    stop("There is no input for prediction after NA omit")

  }else{

  prediction_result <- imp_prediction(final_df)

  }
  return(prediction_result)

}

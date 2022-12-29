#' Mutation impact Prediction
#'
#' Main function for mutation impact prediction
#'
#' With an input data.frame which contains 'PDB_ID', 'Chain', 'Position',
#' 'Reference Amino Acid', 'Mutated Amino Acid' and 'Gene_Name' (optional) information respectively,
#' make a prediction about the impact of a mutation and classifies the mutation as 'Pathogenic' or 'Neutral'.
#' Predictions can be done by using 7 angstrom cutoff - all atom approach or 5 angstrom cutoff - ca atom approach.
#' Other cutoffs can be used for exploratory purposes.
#'
#' @param info_df data.frame containing 'PDB_ID', 'Chain', 'Position',
#' 'Reference Amino Acid', 'Mutated Amino Acid' and 'Gene_Name' (optional) information respectively.
#' @param PDB_path PDB file path (default = NULL)
#' @param n_threads number of threads (default = NULL)
#' @param gene_name_info whether there is gene name information in the input or not (default = TRUE)
#' @param distance_cutoff distance cutoff for setting edges (default = 5)
#' @param network_approach network building approach; "all" (default) for using all atoms only or "ca" for using ca atoms only
#'
#' @return data.frame which contains prediction results
#'
#' @export
#'

predatoR <- function(info_df, PDB_path = NULL, n_threads = NULL, gene_name_info = TRUE, distance_cutoff = 5, network_approach = "all"){

  if(!network_approach %in% c("all", "ca")){
    stop("Network approach needs to be 'all' or 'ca'")

  }

  if(is.null(n_threads) == TRUE){
    n.cores <- parallel::detectCores() - 1
  }else{
    n.cores <- n_threads
  }

  if(is.data.frame(info_df) == F){
    stop("Input should be a data.frame.")
  }
  if(gene_name_info == TRUE){
    if(ncol(info_df) != 6){
      stop("Input data.frame should contain 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name respectively.")
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
      info_df$Gene_Name <- ""
      info_df$Orig_AA <- toupper(info_df$Orig_AA)
      info_df$Mut_AA <- toupper(info_df$Mut_AA)
    }
  }
  if(length(setdiff(unique(unlist(c(info_df[4], info_df[5]))), colnames(blosum_data))) > 0){

    false_name <- setdiff(unique(unlist(c(info_df[4], info_df[5]))), colnames(blosum_data))

    stop(paste0(paste0(false_name, collapse = ", "),
                " couldn't find in the amino acid names (Amino acid names should be 3 letter codes)"))
  }

  my.cluster <- parallel::makeCluster(n.cores, strategy = "sequential")
  doParallel::registerDoParallel(my.cluster)

  final_df <- data.frame()
  for(i in unique(info_df$PDB_ID)){

    filtered_info_df <- info_df[info_df$PDB_ID == i,]

    atom_matrix <- read_PDB(i, PDB_path = PDB_path, network_approach = network_approach)

    if(is.data.frame(atom_matrix) == FALSE){
      next
    }


    message(crayon::white(paste0("PDB ID:", "\t\t\t\t", i, "\n",
                                 "Position(s):", "\t\t\t",
                                 paste0(sapply(split(filtered_info_df$Position, ceiling(seq_along(filtered_info_df$Position)/20)),
                                               function(x) paste(x, collapse = ", ")), collapse = "\n\t\t\t\t"))))


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

    edge_list <- PDB2connections(atom_matrix, filtered_info_df, n_threads = n.cores, single_run = FALSE, distance_cutoff = distance_cutoff)

    filtered_info_df$degree_z_score <- degree_score(edge_list, filtered_info_df)

    filtered_info_df$eigen_z_score <- eigen_centrality_score(edge_list, filtered_info_df)

    filtered_info_df$shortest_path_z <- shorteset_path_score(edge_list, filtered_info_df)

    filtered_info_df$betweenness_scores_z <- betweenness_score(edge_list, filtered_info_df)

    filtered_info_df$clique_z_score <- clique_score(edge_list, filtered_info_df, n_threads = n.cores, single_run = FALSE)

    filtered_info_df$pagerank_z_score <- pagerank_score(edge_list, filtered_info_df)

    filtered_info_df <- gnomad_scores(filtered_info_df)

    filtered_info_df$blosum62_scores <- as.numeric(BLOSUM62_score(filtered_info_df))

    filtered_info_df <- KEGG_pathway_number(filtered_info_df)

    filtered_info_df <- genic_intolerance(filtered_info_df)

    filtered_info_df <- GO_terms(filtered_info_df)

    filtered_info_df <- DisGeNET(filtered_info_df)

    filtered_info_df <- gene_essentiality(filtered_info_df)

    filtered_info_df <- GTEx(filtered_info_df)

    filtered_info_df <- amino_acid_features(filtered_info_df)

    final_df <- rbind(final_df, filtered_info_df)
  }


  if(network_approach == "ca" & distance_cutoff == 7){
    final_df <- stats::na.omit(final_df)

      if(nrow(final_df) == 0){

        parallel::stopCluster(my.cluster)

        stop("There is no input for prediction")

      }else{

        prediction_result <- impact_prediction(final_df, distance_cutoff = distance_cutoff, network_approach = network_approach)

        parallel::stopCluster(my.cluster)

        return(prediction_result)

      }

  }

  else if(network_approach == "all" & distance_cutoff == 5){
    final_df <- stats::na.omit(final_df)

    if(nrow(final_df) == 0){

      parallel::stopCluster(my.cluster)

      stop("There is no input for prediction")

    }else{

      prediction_result <- impact_prediction(final_df, distance_cutoff = distance_cutoff, network_approach = network_approach)

      parallel::stopCluster(my.cluster)

      return(prediction_result)

    }

  }else{
      return(final_df)
    }

}

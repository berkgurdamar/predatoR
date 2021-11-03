#' Genic Intolerance Score
#'
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return Genic Intolerance scores of input genes
#'
#' @export
#'

genic_intolerance <- function(filtered_info_df){

  if(ncol(filtered_info_df) < 6){
    stop("Input data.frame should contain at least 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name ... respectively.")
  }
  else{
    colnames(filtered_info_df)[1:6] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA", "Gene_Name")
  }

  no_name <- c()
  for(i in 1:nrow(filtered_info_df)){

    genic_int_result <- genic_intolerance_data$ALL_0.1.[which(genic_intolerance_data$GENE == filtered_info_df$Gene_Name[i])]

    if(length(genic_int_result) == 0){

      no_name <- c(no_name, filtered_info_df$Gene_Name[i])

      filtered_info_df$genic_intolerance[i] <- NA
    }
    else{

      filtered_info_df$genic_intolerance[i] <- genic_int_result

    }

  }

  if(length(unique(no_name)) > 0){
    message(crayon::white(paste0("\n", "Genic Intolerance score of the gene(s) ", paste0(unique(no_name), collapse = ", "),
                                 " couldn't find, will be removed from the query", "\n")))
  }

  message(crayon::white(paste0("Genic Intolerance Score:", "\t", "DONE")))

  return(filtered_info_df)
}

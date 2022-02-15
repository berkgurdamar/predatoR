#' GO Term Number from Gene Name
#'
#' Function finds the number of GO Terms which related with the input genes.
#'
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return Number of GO Terms which related with the input genes
#'
#' @export
#'

GO_terms <- function(filtered_info_df){
  if(ncol(filtered_info_df) < 6){
    stop("Input data.frame should contain at least 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name ... respectively.")
  }
  else{
    colnames(filtered_info_df)[1:6] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA", "Gene_Name")
  }

  for(i in 1:nrow(filtered_info_df)){
    filtered_info_df$go_terms[i] <- length(go_terms$GO.term.name[which(go_terms$Gene.name == filtered_info_df$Gene_Name[i])])
  }

  message(crayon::white(paste0("GO Term Number:", "\t\t\t", "DONE")))

  return(filtered_info_df)
}

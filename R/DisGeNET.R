#' Disease-Gene relationship from DisGeNET
#'
#' Function finds the number of diseases related with the input genes.
#'
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return Number of diseases related with input gene
#' @export
#'


DisGeNET <- function(filtered_info_df){
  if(ncol(filtered_info_df) < 6){
    stop("Input data.frame should contain at least 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name ... respectively.")
  }
  else{
    colnames(filtered_info_df)[1:6] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA", "Gene_Name")
  }

  for(i in 1:nrow(filtered_info_df)){

    filtered_info_df$disgenet[i] <- length(which(filtered_info_df$Gene_Name[i] == disgenet$geneSymbol))
  }

  message(crayon::white(paste0("DisGeNET Disease Number:", "\t", "DONE")))

  return(filtered_info_df)
}

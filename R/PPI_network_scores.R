#' PPI Network Scores
#'
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return Network property z-scores of input gene name
#' @export
#'

PPI_network_scores <- function(filtered_info_df){

  if(ncol(filtered_info_df) < 6){
    stop("Input data.frame should contain at least 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name ... respectively.")
  }
  else{
    colnames(filtered_info_df)[1:6] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA", "Gene_Name")
  }

  for(i in 1:nrow(filtered_info_df)){
    idx <- which(filtered_info_df$Gene_Name[i] == rownames(PPI_properties))
    if(length(idx) > 0){
      for(j in colnames(PPI_properties))
      filtered_info_df[i,j] <- PPI_properties[idx, j]
    }
  }


  message(crayon::white(paste0("PPI Network Scores:", "\t\t", "DONE")))

  return(filtered_info_df)
}

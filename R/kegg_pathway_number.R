#' KEGG Pathway Number from Gene Name
#'
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return Number of KEGG pathways which contain the input genes
#'
#' @export
#'

KEGG_pathway_number <- function(filtered_info_df){

  if(ncol(filtered_info_df) < 6){
    stop("Input data.frame should contain at least 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name ... respectively.")
  }
  else{
    colnames(filtered_info_df)[1:6] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA", "Gene_Name")
  }

  for(i in 1:nrow(filtered_info_df)){
    filtered_info_df$kegg_pathway_number[i] <- length(grep(filtered_info_df$Gene_Name[i], kegg_info))
  }

  message(crayon::white(paste0("KEGG Pathway Number:", "\t\t", "DONE")))

  return(filtered_info_df)

}

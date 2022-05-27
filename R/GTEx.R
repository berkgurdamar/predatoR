#' GTEx Median Expression Function
#'
#' `GTEx()` function calculates the median gene expression value with using 54 different tissue types' median gene expression value from
#' [GTEx](https://gtexportal.org/home/) and returns the value of input gene.
#'
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return median gene expression values of input genes
#' @export
#'


GTEx <- function(filtered_info_df){
  if(ncol(filtered_info_df) < 6){
    stop("Input data.frame should contain at least 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name ... respectively.")
  }
  else{
    colnames(filtered_info_df)[1:6] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA", "Gene_Name")
  }

  gtex_median <- gtex_median[which(toupper(gtex_median$gene_name) %in% toupper(filtered_info_df$Gene_Name)),]

  for(i in 1:nrow(filtered_info_df)){

    idx <- which(toupper(gtex_median$gene_name) %in% toupper(filtered_info_df$Gene_Name[i]))

    if(length(idx) > 0){
      filtered_info_df$gtex[i] <- max(gtex_median$median[idx])
    }else{
      filtered_info_df$gtex[i] <- NA
    }
  }

  message(crayon::white(paste0("GTEx Score:", "\t\t\t", "DONE")))

  return(filtered_info_df)
}

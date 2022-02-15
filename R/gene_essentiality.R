#' Gene Essentiality Score
#'
#' Function finds the gene essentiality scores of the input genes from OGEE database.
#'
#' @param filtered_info_df input data.frame which contain only one PDB entries
#'
#' @return Gene Essentiality Scores of input genes
#' @export
#'



gene_essentiality <- function(filtered_info_df){
  if(ncol(filtered_info_df) < 6){
    stop("Input data.frame should contain at least 6 columns; PDB_ID, Chain, Position, Orig_AA, Mut_AA and Gene_Name ... respectively.")
  }
  else{
    colnames(filtered_info_df)[1:6] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA", "Gene_Name")
  }

  gene_essentiality_scores <- gene_essentiality_scores[which(toupper(gene_essentiality_scores$gene) %in% toupper(filtered_info_df$Gene_Name)),]

  for(i in 1:nrow(filtered_info_df)){

    idx <- which(toupper(gene_essentiality_scores$gene) %in% toupper(filtered_info_df$Gene_Name[i]))

    if(length(idx) > 0){
      filtered_info_df$gene_essentiality[i] <- max(gene_essentiality_scores$score[idx])
    }else{
      filtered_info_df$gene_essentiality[i] <- NA
    }
  }

  message(crayon::white(paste0("Gene Essentiality Score:", "\t", "DONE")))

  return(filtered_info_df)
}

#' BLOSUM62 Score
#'
#' @param info_df data.frame contains all the input mutations
#'
#' @return BLOSUM62 score of input mutation
#'
#' @export
#'

BLOSUM62_score <- function(info_df){

  colnames(info_df)[1:5] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA")

  blosum62_scores <- c()
  for(i in 1:nrow(info_df)){
    blosum62_scores <- c(blosum62_scores, blosum_data[info_df$Orig_AA[i], info_df$Mut_AA[i]])
  }
  message(crayon::white(paste0("BLOSUM62 Score:", "\t\t\t", "DONE")))

  return(blosum62_scores)
}



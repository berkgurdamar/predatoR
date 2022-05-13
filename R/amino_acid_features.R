#' Amino Acid Features
#'
#' This function adds accessible surface area and hydrophobicity scale
#' values of reference amino acid, mutant amino acid and difference between mutant and reference
#' amino acids.
#'
#' @param info_df data.frame contains all the input mutations
#'
#' @return accessible surface area and hydrophobicity scale
#' values of reference amino acid, mutant amino acid and difference between mutant and reference
#' amino acids
#' @export
#'

amino_acid_features <- function(info_df){

  colnames(info_df)[1:5] <- c("PDB_ID", "Chain", "Position", "Orig_AA", "Mut_AA")

  for(i in 1:nrow(info_df)){
    info_df$ref_asa[i] <- amino_acid_features_df[1,info_df$Orig_AA[i]]
    info_df$mut_asa[i] <- amino_acid_features_df[1,info_df$Mut_AA[i]]
    info_df$asa_diff[i] <- amino_acid_features_df[1,info_df$Mut_AA[i]] - amino_acid_features_df[1,info_df$Orig_AA[i]]

    info_df$ref_hyd[i] <- amino_acid_features_df[2,info_df$Orig_AA[i]]
    info_df$mut_hyd[i] <- amino_acid_features_df[2,info_df$Mut_AA[i]]
    info_df$hyd_diff[i] <- amino_acid_features_df[2,info_df$Mut_AA[i]] - amino_acid_features_df[2,info_df$Orig_AA[i]]

  }

  message(crayon::white(paste0("Amino Acid Features:", "\t\t", "DONE")))

  return(info_df)
}

#' Mutation Type Classifier
#'
#' This functios classifies the Charged-Charged, Charged-Polar, Polar-Charged, Polar-Polar and
#' Hydrophobic-Hydrophobic mutations as 'Not deadly (0)', classifies Charged-Hydrophobic,
#' Hydrophobic-Charged, Hydrophobic-Polar, Polar-Hydrophobic mutations as 'Deadly (1)'.
#'
#' @param info_df data.frame contains all the input mutations
#'
#' @return Mutation Type as a factor; 1 = Deadly, 0 = Not deadly
#' @export
#'

mutation_type <- function(info_df){

  mut_type <- c()
  for(i in 1:nrow(info_df)){

    amino_acids <- toupper(c("ala", "arg", "asn", "asp", "val",
                             "cys", "glu", "gln", "gly", "tyr",
                             "his", "ile", "leu", "lys", "met",
                             "phe", "pro", "ser", "thr", "trp"))

    aa_properties <- c("H", "C", "P", "C", "H",
                       "P", "C", "P", "P", "P",
                       "C", "H", "H", "C", "H",
                       "H", "H", "P", "P", "H")


    mut_type <- c(mut_type, paste0(aa_properties[which(info_df$Orig_AA[i] == amino_acids)],
                                   aa_properties[which(info_df$Mut_AA[i] == amino_acids)]))

  }
  mut_type[which(mut_type %in% c("CH", "HC", "HP", "PH"))] <- 1
  mut_type[which(mut_type %in% c("CC",  "CP", "HH", "PC", "PP"))] <- 0

  message(crayon::white(paste0("Mutation Type:", "\t\t\t", "DONE")))

  return(mut_type)
}

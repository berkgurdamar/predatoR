#' Genic Intolerance Score
#'
#' @param gene_name Gene name
#'
#' @return numeric
#' @export
#'

genic_intolerance <- function(gene_name){

  genic_int_result <- genic_intolerance_data$ALL_0.1.[which(genic_intolerance_data$GENE == gene_name)]

  if(is.na(genic_int_result) == T){
    message(crayon::white(paste0("\n", "Genic Intolerance score of the gene ", gene_name,  " couldn't find, will be removed from the query", "\n")))

    message(crayon::white(paste0("Genic Intolerance Score:", "\t", "DONE")))

    return(NA)
  }

  message(crayon::white(paste0("Genic Intolerance Score:", "\t", "DONE")))

  return(genic_int_result)

}

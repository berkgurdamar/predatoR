#' KEGG Pathway Number from Gene Name
#'
#' @param gene_name Gene name
#'
#' @return Number of KEGG pathways which contain the input gene
#' @export
#'

KEGG_pathway_number <- function(gene_name){

  message(crayon::white(paste0("Number of KEGG Pathways:", "\t", "DONE")))

  return(length(grep(gene_name, kegg_info)))

}

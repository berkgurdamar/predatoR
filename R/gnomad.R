#' GnomAD Scores
#'
#' @param PDB_ID PDB ID
#' @param info_df data.frame contains all the input mutations
#'
#' @return list that contains a data.frame containing GnomAD scores of input mutations and
#' the gene name
#'
#' @export
#'

gnomad <- function(PDB_ID, info_df){


  res <- biomart_data[which(biomart_data$pdb == PDB_ID),]

  if(length(res$external_gene_name) > 1){
    gene_name <- readline(prompt = paste0("There are multiple gene names for PDB ", PDB_ID,
                                          ", please write the input gene name (",
                                          paste0(res$external_gene_name, collapse = ", "),
                                          "): "))

    if(!(gene_name %in% res$external_gene_name)){

      stop(paste0("Input gene name ", gene_name,
                  " not found in query genes (",
                  paste0(res$external_gene_name, collapse = ", "), ")"))
    }

    idx <- which(gnomad_data$gene == gene_name)

    if(length(idx) == 0){
      message(crayon::white(paste0("\n", "gnomAD scores of the gene ", gene_name,
                                   " couldn't find, will be removed from the query", "\n")))

      info_df$syn_z <- NA

      info_df$mis_z <- NA

      info_df$pLI <- NA

      message(crayon::white(paste0("GNOMAD Information:", "\t\t", "DONE")))

      return(list(info_df, gene_name))
    }
    else{

      info_df$syn_z <- as.numeric(gnomad_data$syn_z[idx])

      info_df$mis_z <- as.numeric(gnomad_data$mis_z[idx])

      info_df$pLI <- as.numeric(gnomad_data$pLI[idx])

      message(crayon::white(paste0("GNOMAD Information:", "\t\t", "DONE")))

      return(list(info_df, gene_name))
    }

  }
  else if(length(res$external_gene_name) == 0){
    gene_name <- readline(prompt = paste0("No gene name found for the PDB ", PDB_ID,
                                          ", please write the input gene name: "))

    idx <- which(gnomad_data$gene == gene_name)

    if(length(idx) == 0){
      message(crayon::white(paste0("\n", "gnomAD scores of the gene ", gene_name,
                                   " couldn't find, will be removed from the query", "\n")))

      info_df$syn_z <- NA

      info_df$mis_z <- NA

      info_df$pLI <- NA

      message(crayon::white(paste0("GNOMAD Information:", "\t\t", "DONE")))

      return(list(info_df, gene_name))
    }
    else{

      info_df$syn_z <- as.numeric(gnomad_data$syn_z[idx])

      info_df$mis_z <- as.numeric(gnomad_data$mis_z[idx])

      info_df$pLI <- as.numeric(gnomad_data$pLI[idx])

      message(crayon::white(paste0("GNOMAD Information:", "\t\t", "DONE")))

      return(list(info_df, gene_name))
    }
  }
  else{
    idx <- which(gnomad_data$gene == res$external_gene_name)

    info_df$syn_z <- as.numeric(gnomad_data$syn_z[idx])

    info_df$mis_z <- as.numeric(gnomad_data$mis_z[idx])

    info_df$pLI <- as.numeric(gnomad_data$pLI[idx])

    message(crayon::white(paste0("GNOMAD Information:", "\t\t", "DONE")))

    return(list(info_df, gnomad_data$gene[idx]))
  }
}


